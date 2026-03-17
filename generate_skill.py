import asyncio
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import tyro
from anthropic import AsyncAnthropic
from dotenv import load_dotenv
from upskill.generate import GENERATION_PROMPT, generate_skill, parse_json_response
from upskill.models import Skill, SkillMetadata

from src import constants, utils
from src.constants import SKILL_MODEL
from src.prompts.default import TASK_DESCRIPTION
from src.agents import MDAgent, PrepAgent
from src.upskill.md_evaluator import MDEvalResult, evaluate_md_run, parse_system_name
from src.upskill.trace_exporter import export_run, extract_error_examples

load_dotenv(dotenv_path=constants.ENV_FILE)

# Module-level logger — initially stdout-only; main() re-points it to the
# appropriate log file once the output directory is known.
logger = utils.get_class_logger("GenerateSkill", log_to_file=False)


def _setup_logger(log_file: Path) -> None:
    """Re-configure the module logger to also write to *log_file*."""
    global logger
    # Drop existing handlers to avoid duplicates when called multiple times.
    logger.handlers.clear()
    logger = utils.get_class_logger("GenerateSkill", log_file=log_file)

# System prompt for MD-specific refinement — removes the 200-400 word limit
# from upskill's default GENERATION_PROMPT since MD skills require more detail.
_MD_REFINEMENT_PROMPT = (
    GENERATION_PROMPT.replace(
        "body: 200-400 word markdown guide with step-by-step instructions and 2-3 examples",
        "body: detailed markdown guide covering the full pipeline, specific error-recovery "
        "patterns, and actionable fixes for the reported failures — be thorough",
    )
)


@dataclass
class Args:
    run_index: int = -1
    "Which run to export from the agent log (0 = first, -1 = last)."

    skill_name: str = ""
    "Name for the generated skill directory. Defaults to 'md-run-<timestamp>'."

    eval_model: list[str] = field(default_factory=list)
    "Model(s) to run legacy upskill Q&A eval on after generation (e.g. haiku, sonnet)."

    eval_only: bool = False
    "Skip skill generation; only run evaluation on an existing skill."

    eval_runs: int = 3
    "Number of evaluation runs per model (for the legacy Q&A eval)."

    trace_output: Path = constants.AGENT_LOGS / "agent-trace.txt"
    "Where to write the exported trace file."

    eval_systems: list[str] = field(default_factory=list)
    "MD pipeline runs to score with file-based metrics. Format: 'PDBID_LIGAND:path/to/sandbox'."

    compare_systems: list[str] = field(default_factory=list)
    "Systems to run baseline-vs-skilled comparison on. Format: 'PDBID_LIGAND' or 'PDBID_None'."

    compare_model: str = constants.MODEL_NAME
    "Model to use for compare_systems runs."

    compare_temp: float | None = None
    "Simulation temperature in Kelvin for compare_systems runs."

    compare_duration: float | None = None
    "Simulation duration in nanoseconds for compare_systems runs."


async def _generate_skill_async(examples: list[str]) -> Skill:
    """Generate a skill from the task description and error-recovery examples."""
    return await generate_skill(
        task=TASK_DESCRIPTION,
        examples=examples or None,
        model=SKILL_MODEL,
    )


async def _refine_md_skill_async(skill: Skill, failures: list[str]) -> Skill:
    """
    Refine an MD skill based on pipeline step failures.

    Uses a custom prompt that passes the full skill body (upskill's built-in
    refine_skill() truncates to 500 chars, which loses most of an MD skill).
    """
    client = AsyncAnthropic()

    prompt = (
        f"Improve this MD simulation skill based on pipeline execution failures.\n\n"
        f"Name: {skill.name}\n"
        f"Description: {skill.description}\n\n"
        f"Current skill body:\n{skill.body}\n\n"
        f"Pipeline step failures (steps that did not produce expected output files):\n"
        + "\n".join(f"- {f}" for f in failures)
        + "\n\nOutput the improved skill as JSON (same structure, no code blocks)."
    )

    response = await client.messages.create(
        model=SKILL_MODEL,
        max_tokens=8192,
        system=_MD_REFINEMENT_PROMPT,
        messages=[{"role": "user", "content": prompt}],
    )

    data = parse_json_response(response.content[0].text)
    return Skill(
        name=data.get("name", skill.name),
        description=data.get("description", skill.description),
        body=data["body"],
        references=data.get("references", skill.references),
        scripts=data.get("scripts", skill.scripts),
        metadata=SkillMetadata(
            generated_by=SKILL_MODEL,
            generated_at=datetime.now(timezone.utc),
            source_task=skill.metadata.source_task,
        ),
    )


def _load_skill_from_dir(skill_dir: Path) -> Skill | None:
    """Load an existing Skill object from a skill directory."""
    from src.upskill.skill_loader import load_skill_object
    return load_skill_object(skill_dir)


def _md_eval_to_failures(result: MDEvalResult) -> list[str]:
    """Convert an MDEvalResult into human-readable failure strings for refine."""
    return [
        f"Pipeline step '{step}' incomplete (score {score:.2f}) — "
        f"expected output files were missing or empty"
        for step, score in result.step_scores.items()
        if score < 1.0
    ]


def _run_agent_and_score(
    system_name: str,
    model_name: str,
    use_skills: bool,
    sandbox_dir: Path,
    md_temp: float | None = None,
    md_duration: float | None = None,
) -> MDEvalResult:
    """Run the full PrepAgent + MDAgent pipeline and score output files."""
    pdb_id, ligand_name = parse_system_name(system_name)

    sandbox_dir.mkdir(parents=True, exist_ok=True)
    log_file = sandbox_dir / "run.log"
    label = "skilled" if use_skills else "baseline"
    logger = utils.get_class_logger(f"CompareEval({label})", log_file=log_file)
    logger.info(f"Starting {label} run for {system_name}")

    prep_agent = PrepAgent(
        model_name=model_name,
        temperature=constants.TEMPERATURE,
        sandbox_dir=sandbox_dir,
        pdb_id=pdb_id,
        ligand_name=ligand_name,
        md_temp=md_temp,
        md_duration=md_duration,
    )
    prep_agent.setup_tools()
    pdb_file_path, resolved_ligand, plan, _ = prep_agent.run()

    md_agent = MDAgent(
        model_name=model_name,
        temperature=constants.TEMPERATURE,
        sandbox_dir=sandbox_dir,
        structure_path=pdb_file_path,
        pdb_id=pdb_id,
        ligand_name=resolved_ligand,
        md_temp=plan["parameters"]["temperature_k"],
        md_duration=plan["parameters"]["duration_ns"],
        plan=plan,
        use_skills=use_skills,
        post_pipeline_chat=False,
    )
    md_agent.setup_tools()
    md_agent.run()

    return evaluate_md_run(sandbox_dir, system_name)


def _log_comparison(system_name: str, baseline: MDEvalResult, skilled: MDEvalResult) -> None:
    system_type = "protein-ligand" if skilled.is_protein_ligand else "protein-only"
    all_steps = list(skilled.step_scores.keys())
    col = 32

    lines = [
        f"System: {system_name} ({system_type})",
        f"  {'Step':<{col}} {'BASELINE':>10}  {'WITH SKILL':>10}  {'DELTA':>8}",
        f"  {'-' * (col + 34)}",
    ]
    for step in all_steps:
        b = baseline.step_scores.get(step, 0.0)
        s = skilled.step_scores.get(step, 0.0)
        delta = s - b
        sign = "+" if delta > 0 else ""
        lines.append(f"  {step:<{col}} {b:>10.2f}  {s:>10.2f}  {sign}{delta:>7.2f}")
    lines.append(f"  {'-' * (col + 34)}")
    b_total = baseline.overall_score
    s_total = skilled.overall_score
    delta_total = s_total - b_total
    sign = "+" if delta_total > 0 else ""
    lines.append(f"  {'Overall':<{col}} {b_total:>10.2f}  {s_total:>10.2f}  {sign}{delta_total:>7.2f}")
    lines.append(f"  {'Pipeline success':<{col}} {str(baseline.pipeline_successful):>10}  {str(skilled.pipeline_successful):>10}")
    logger.info("\n".join(lines))


def _log_md_eval_result(result: MDEvalResult) -> None:
    system_type = "protein-ligand" if result.is_protein_ligand else "protein-only"
    lines = [f"System: {result.system_name} ({system_type})"]
    for step, score in result.step_scores.items():
        lines.append(f"  {step:<30}: {score:.2f}")
    lines.append(f"  {'Overall':<30}: {result.overall_score:.2f} | SUCCESS: {result.pipeline_successful}")
    logger.info("\n".join(lines))


def main(args: Args) -> None:
    skill_name = args.skill_name or f"md-run-{utils.time_now()}"
    skill_dir = constants.SKILLS_DIR / skill_name

    if not args.eval_only:
        # Log into the skill directory alongside SKILL.md
        skill_dir.mkdir(parents=True, exist_ok=True)
        _setup_logger(skill_dir / "generate.log")

        # Step 1: export trace (kept for reference / manual inspection)
        logger.info(f"Exporting run {args.run_index} from {constants.JSON_LOG_FILE} ...")
        trace_path = export_run(run_index=args.run_index, output_path=args.trace_output)
        logger.info(f"Trace written to: {trace_path}")

        # Step 2: extract error-recovery examples and generate skill via Python API
        error_examples = extract_error_examples(run_index=args.run_index)
        logger.info(f"Extracted {len(error_examples)} error-recovery examples from trace.")
        if error_examples:
            for ex in error_examples:
                logger.info(f"  • {ex[:100]}{'...' if len(ex) > 100 else ''}")

        logger.info(f"Generating skill '{skill_name}' with {SKILL_MODEL} ...")
        try:
            skill = asyncio.run(_generate_skill_async(error_examples))
        except Exception as e:
            logger.error(f"Skill generation failed: {e}")
            sys.exit(1)

        skill.save(skill_dir)
        skill_md = skill_dir / "SKILL.md"
        tokens_approx = int(len(skill.body.split()) * 1.3)
        logger.info(f"Skill saved to: {skill_md}  (~{tokens_approx} tokens)")
        all_skills = sorted(constants.SKILLS_DIR.iterdir()) if constants.SKILLS_DIR.exists() else []
        logger.info(f"Total skills in library: {len([d for d in all_skills if (d / 'SKILL.md').exists()])}")

    # Step 3 (optional): legacy upskill Q&A eval — note this tests declarative
    # knowledge, not pipeline execution. Use --compare-systems for real evaluation.
    for model in args.eval_model:
        logger.info(f"Evaluating skill on model: {model} (legacy Q&A eval)")
        import subprocess
        rc = subprocess.run([
            "upskill", "eval", str(skill_dir),
            "--model", model,
            "--runs", str(args.eval_runs),
        ]).returncode
        if rc != 0:
            logger.warning(f"upskill eval exited with code {rc} for model {model}.")

    # Step 4 (optional): score an existing run's output files without re-running
    if args.eval_systems:
        logger.info("MD pipeline file-based evaluation")
        for entry in args.eval_systems:
            if ":" not in entry:
                logger.error(f"Invalid --eval-systems entry '{entry}'. Expected 'SYSTEM_NAME:path'.")
                continue
            system_name, sandbox_path = entry.split(":", 1)
            result = evaluate_md_run(sandbox_path, system_name)
            _log_md_eval_result(result)

    # Step 5 (optional): baseline vs. skilled comparison, with automatic refinement
    # if the skilled run does not complete the pipeline successfully.
    if args.compare_systems:
        logger.info("MD baseline vs skilled comparison")
        for system_name in args.compare_systems:
            run_id = utils.time_now()
            base_dir = constants.DATA_DIR / f"compare_{run_id}_{system_name}"

            # Log comparison orchestration into the comparison base directory
            base_dir.mkdir(parents=True, exist_ok=True)
            _setup_logger(base_dir / "eval.log")

            logger.info(f"[{system_name}] Running BASELINE (no skill) ...")
            baseline = _run_agent_and_score(
                system_name=system_name,
                model_name=args.compare_model,
                use_skills=False,
                sandbox_dir=base_dir / "baseline",
                md_temp=args.compare_temp,
                md_duration=args.compare_duration,
            )

            logger.info(f"[{system_name}] Running WITH SKILL ...")
            skilled = _run_agent_and_score(
                system_name=system_name,
                model_name=args.compare_model,
                use_skills=True,
                sandbox_dir=base_dir / "skilled",
                md_temp=args.compare_temp,
                md_duration=args.compare_duration,
            )

            _log_comparison(system_name, baseline, skilled)

            # Refine skill if the pipeline did not complete successfully.
            if not skilled.pipeline_successful:
                failures = _md_eval_to_failures(skilled)
                if failures:
                    logger.info(f"[{system_name}] Skill run incomplete — refining skill based on "
                                f"{len(failures)} failed step(s) ...")
                    existing_skill = _load_skill_from_dir(skill_dir)
                    if existing_skill:
                        try:
                            refined = asyncio.run(_refine_md_skill_async(existing_skill, failures))
                            refined.save(skill_dir)
                            tokens_approx = int(len(refined.body.split()) * 1.3)
                            logger.info(f"Refined skill saved to {skill_dir / 'SKILL.md'} "
                                        f"(~{tokens_approx} tokens)")
                        except Exception as e:
                            logger.warning(f"Skill refinement failed: {e}")
                    else:
                        logger.warning(f"Could not load skill from {skill_dir} for refinement.")
            else:
                logger.info(f"[{system_name}] Skilled pipeline successful — no refinement needed.")


if __name__ == "__main__":
    args = tyro.cli(Args)
    main(args)

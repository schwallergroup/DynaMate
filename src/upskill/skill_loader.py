from pathlib import Path
from upskill.models import Skill
from src import constants


def load_skill_object(skill_dir: Path) -> Skill | None:
    """
    Parse a skill directory into an upskill Skill object.

    Handles both the current upskill render() format (# header) and the
    older YAML-frontmatter format produced by earlier versions.

    Returns None if no SKILL.md exists or parsing fails.
    """
    skill_md = skill_dir / "SKILL.md"
    if not skill_md.exists():
        return None

    content = skill_md.read_text(encoding="utf-8").strip()
    if not content:
        return None

    try:
        lines = content.split("\n")
        name = lines[0].lstrip("# ").strip() if lines else skill_dir.name
        description = lines[2].strip() if len(lines) > 2 else ""
        # Skip the "## Instructions" header (and trailing blank line) that
        # render() adds, so it won't be duplicated on the next save().
        body_start = 4
        if len(lines) > body_start and lines[body_start].strip().lower() == "## instructions":
            body_start += 1
            if len(lines) > body_start and lines[body_start].strip() == "":
                body_start += 1
        body = "\n".join(lines[body_start:]).strip() if len(lines) > body_start else content

        # Sanitize name to meet upskill's validation requirement.
        import re
        name = re.sub(r"[^a-z0-9-]", "-", name.lower())
        name = re.sub(r"-+", "-", name).strip("-") or "skill"

        return Skill(name=name, description=description, body=body)

    except Exception:
        return None


def load_all_skills() -> str:
    """
    Load every skill found under SKILLS_DIR and concatenate their bodies
    for injection into agent system prompts.

    Skills are loaded newest-first (by directory mtime) so the most
    recently generated examples appear first. Total injected content is
    capped at MAX_SKILL_CHARS to avoid blowing the context window.
    """
    if not constants.SKILLS_DIR.exists():
        return ""

    skill_dirs = sorted(
        [d for d in constants.SKILLS_DIR.iterdir() if (d / "SKILL.md").exists()],
        key=lambda d: d.stat().st_mtime,
        reverse=True,
    )

    parts = []
    total_chars = 0

    for skill_dir in skill_dirs:
        skill = load_skill_object(skill_dir)
        if not skill or not skill.body:
            continue
        block = f'<skill name="{skill.name}">\n{skill.body}\n</skill>\n\n'
        if total_chars + len(block) > constants.MAX_SKILL_CHARS:
            break
        parts.append(block)
        total_chars += len(block)

    return "".join(parts)

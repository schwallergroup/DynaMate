import re
import csv
from datetime import datetime
from pathlib import Path

LOG_DIR = Path("agent_logs")
OUT_FILE = "run_summary.csv"
MODEL_NAME = "claude-3-5-haiku-20241022"  # update this for each model run


def parse_timestamp(line):
    try:
        return datetime.strptime(line.split(" - ")[0], "%Y-%m-%d %H:%M:%S,%f")
    except Exception:
        return None


def parse_logs():
    prep_log = (LOG_DIR / "PrepAgent.log").read_text()
    md_log = (LOG_DIR / "MDAgent.log").read_text()

    # --- Step 1: Identify runs and their metadata ---
    runs = []
    current_run = {}
    for line in prep_log.splitlines():
        if "PrepAgent initialized." in line:
            if current_run:
                runs.append(current_run)
            current_run = {"tools": set(), "start_time": parse_timestamp(line)}
        elif "User input:" in line:
            current_run["protein"] = line.split("User input:")[1].strip()
        elif "User requested ligand:" in line:
            current_run["ligand"] = line.split("User requested ligand:")[1].strip()
        elif "Executing tool:" in line:
            tool_match = re.search(r"Executing tool:\s*(\w+)", line)
            if tool_match:
                current_run["tools"].add(tool_match.group(1))
    if current_run:
        runs.append(current_run)

    # --- Step 2: Parse MDAgent log ---
    md_lines = md_log.splitlines()
    md_runs = []
    current = None
    for line in md_lines:
        if "MDAgent initialized." in line:
            if current:
                md_runs.append(current)
            current = {
                "start_time": parse_timestamp(line),
                "iterations": 0,
                "tools_called": set(),
                "attempted": 0,
                "success": 0,
            }
        elif "Logging agent iteration" in line and current:
            current["iterations"] += 1
        elif "Executing tool:" in line and current:
            current["attempted"] += 1
            match = re.search(r"Executing tool:\s*(\w+)", line)
            if match:
                current["tools_called"].add(match.group(1))
        elif "Tool result:" in line and current:
            if "Error" not in line:
                current["success"] += 1
        elif "MD Pipeline completed successfully" in line and current:
            current["end_time"] = parse_timestamp(line)
    if current:
        md_runs.append(current)

    # --- Step 3: Merge PrepAgent + MDAgent runs by order ---
    rows = []
    for i, run in enumerate(runs):
        md = md_runs[i] if i < len(md_runs) else {}
        start = md.get("start_time")
        end = md.get("end_time")
        total_time = (end - start).total_seconds() if start and end else None
        rows.append(
            {
                "Model": MODEL_NAME,
                "Protein": run.get("protein", "Unknown"),
                "Ligand": run.get("ligand", "None"),
                "Iterations": md.get("iterations", 0),
                "Total Time (s)": total_time,
                "Subtasks Attempted": md.get("attempted", 0),
                "Subtasks Successful": md.get("success", 0),
                "Tools Called": ", ".join(sorted(md.get("tools_called", run["tools"]))),
            }
        )

    # --- Step 4: Write to CSV ---
    fieldnames = [
        "Model",
        "Protein",
        "Ligand",
        "Iterations",
        "Total Time (s)",
        "Subtasks Attempted",
        "Subtasks Successful",
        "Tools Called",
    ]

    with open(OUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Summary written to {OUT_FILE}")


if __name__ == "__main__":
    parse_logs()

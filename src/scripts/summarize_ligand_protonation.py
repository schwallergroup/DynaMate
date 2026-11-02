"""
Summarizes the ligand protonation records (*_protonation.json) of many DynaMate runs.

Usage:
    python -m src.scripts.summarize_ligand_protonation <runs_dir> [--csv out.csv]

Reports, over all runs: agreement between Dimorphite-DL and pkasolver, ambiguous states, bond-order source,
how often the agent changed the default state, and the manual overrides (states predicted by neither tool)
that should be reviewed by hand.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def collect(runs_dir: Path) -> list[dict]:
    rows = []
    for path in sorted(runs_dir.rglob("*_protonation.json")):
        if "stale_protonation_" in str(path):
            continue
        record = json.loads(path.read_text())
        selected = record["selected"]
        default = record["history"][0] if record["history"] else selected
        rows.append({
            "run": path.parent.name,
            "ligand": record["ligand_name"],
            "ph": record["ph"],
            "bond_order_source": record["bond_order_source"],
            "comparison": record["comparison"],
            "ambiguous": record["ambiguous"],
            "pkasolver_error": record["pkasolver"].get("error") or "",
            "default_smiles": default["smiles"],
            "default_charge": default["net_charge"],
            "selected_smiles": selected["smiles"],
            "selected_charge": selected["net_charge"],
            "selected_source": selected["source"],
            "changed_by_agent": bool(record["history"]),
            "manual_override": record["manual_override"],
            "rationale": selected["rationale"],
            "reference": selected["reference"] or "",
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs_dir", type=Path)
    parser.add_argument("--csv", type=Path, default=None)
    args = parser.parse_args()

    rows = collect(args.runs_dir)
    if not rows:
        print(f"No *_protonation.json found in {args.runs_dir}")
        return

    n = len(rows)
    print(f"{n} ligand protonation records\n")
    for key in ("comparison", "bond_order_source", "selected_source"):
        print(f"{key}:")
        for value, count in Counter(r[key] for r in rows).most_common():
            print(f"  {value:<30} {count:>5}  ({100 * count / n:.1f}%)")
    for key in ("ambiguous", "changed_by_agent", "manual_override"):
        count = sum(r[key] for r in rows)
        print(f"{key:<32} {count:>5}  ({100 * count / n:.1f}%)")
    charge_changed = sum(r["changed_by_agent"] and r["default_charge"] != r["selected_charge"] for r in rows)
    print(f"{'agent changed the net charge':<32} {charge_changed:>5}  ({100 * charge_changed / n:.1f}%)")

    overrides = [r for r in rows if r["manual_override"]]
    if overrides:
        print("\nManual overrides to review:")
        for r in overrides:
            print(f"  {r['run']} {r['ligand']}: {r['default_smiles']} -> {r['selected_smiles']} | {r['rationale']} | {r['reference']}")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nPer-run table written to {args.csv}")


if __name__ == "__main__":
    main()

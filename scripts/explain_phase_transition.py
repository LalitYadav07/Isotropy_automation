from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from _isotropy_paths import OUTPUTS_DIR
from phase_transition_classifier import classify_transition


EXPLAIN_RUNS_DIR = OUTPUTS_DIR / "explain_runs"
EXPLAIN_RUNS_DIR.mkdir(parents=True, exist_ok=True)


def run_command(command: list[str]) -> dict[str, object]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    return {
        "command": command,
        "return_code": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify and route a phase-transition analysis request.")
    parser.add_argument("--parent-cif")
    parser.add_argument("--child-cif")
    parser.add_argument("--distortion-file")
    parser.add_argument("--label")
    parser.add_argument("--include-coupled-catalog", action="store_true")
    parser.add_argument("--max-modeled", type=int, default=60)
    parser.add_argument("--dry-run", action="store_true", help="Classify and write the explain report without executing a downstream workflow.")
    args = parser.parse_args()

    label = args.label or "phase_transition_run"
    run_dir = EXPLAIN_RUNS_DIR / label
    run_dir.mkdir(parents=True, exist_ok=True)

    classification = classify_transition(
        parent_cif=args.parent_cif,
        child_cif=args.child_cif,
        distortion_file=args.distortion_file,
    )

    execution: dict[str, object] | None = None
    if args.dry_run:
        execution = {
            "command": [],
            "return_code": 0,
            "stdout": "",
            "stderr": "",
            "note": "Dry run requested; workflow classification completed without dispatching a downstream script.",
        }
    elif classification.recommended_workflow == "discover_distortion_signatures" and args.parent_cif:
        command = [
            "python3",
            "scripts/discover_distortion_signatures.py",
            args.parent_cif,
            "--label",
            label,
            "--max-modeled",
            str(args.max_modeled),
        ]
        if args.include_coupled_catalog:
            command.append("--include-coupled-catalog")
        execution = run_command(command)
    elif classification.recommended_workflow == "reconstructive_common_subgroups" and args.parent_cif and args.child_cif:
        execution = run_command(
            [
                "python3",
                "scripts/reconstructive_transition_search.py",
                args.parent_cif,
                args.child_cif,
                "--label",
                label,
            ]
        )
    elif classification.recommended_workflow in {
        "commensurate_magnetic_parent_only_discovery",
        "magnetic_parent_child_decomposition",
        "incommensurate_magnetic_superspace_analysis",
        "mixed_magnetic_structural_coupled_analysis",
    }:
        command = [
            "python3",
            "scripts/magnetic_workflows.py",
            "--workflow",
            classification.recommended_workflow,
            "--label",
            label,
        ]
        if args.parent_cif:
            command.extend(["--parent-cif", args.parent_cif])
        if args.child_cif:
            command.extend(["--child-cif", args.child_cif])
        if args.distortion_file:
            command.extend(["--distortion-file", args.distortion_file])
        execution = run_command(command)
    elif classification.recommended_workflow == "generate_tutorial_reports_or_distortion_explainer" and args.distortion_file:
        execution = run_command(
            [
                "python3",
                "scripts/generate_tutorial_reports.py",
                "--input",
                args.distortion_file,
                "--output",
                str(run_dir / "distortion_reports"),
            ]
        )

    report = {
        "label": label,
        "classification": classification.as_dict(),
        "execution": execution,
    }
    (run_dir / "explain_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

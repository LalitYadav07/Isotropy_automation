from __future__ import annotations

import argparse
import json
from pathlib import Path

from _isotropy_paths import OUTPUTS_DIR, resolve_input_path
from transition_metadata import inspect_cif_metadata, inspect_distortion_metadata


MAGNETIC_RUNS_DIR = OUTPUTS_DIR / "magnetic_runs"
MAGNETIC_RUNS_DIR.mkdir(parents=True, exist_ok=True)


WORKFLOW_TITLES = {
    "commensurate_magnetic_parent_only_discovery": "Commensurate magnetic parent-only discovery",
    "magnetic_parent_child_decomposition": "Magnetic parent+child decomposition",
    "incommensurate_magnetic_superspace_analysis": "Incommensurate magnetic superspace analysis",
    "mixed_magnetic_structural_coupled_analysis": "Mixed magnetic-structural coupled analysis",
}


def build_summary(*, workflow: str, label: str, parent_cif: str | None, child_cif: str | None, distortion_file: str | None) -> dict[str, object]:
    parent_path = resolve_input_path(parent_cif) if parent_cif else None
    child_path = resolve_input_path(child_cif) if child_cif else None
    distortion_path = resolve_input_path(distortion_file) if distortion_file else None

    metadata = []
    if parent_path:
        metadata.append(inspect_cif_metadata(parent_path).as_dict())
    if child_path:
        metadata.append(inspect_cif_metadata(child_path).as_dict())
    if distortion_path:
        metadata.append(inspect_distortion_metadata(distortion_path).as_dict())

    notes = {
        "commensurate_magnetic_parent_only_discovery": [
            "Routes magnetic parent-only cases away from structural irrep discovery.",
            "This wrapper is the entry point for future magnetic irrep enumeration and magnetic isotropy subgroup search.",
        ],
        "magnetic_parent_child_decomposition": [
            "Routes magnetic parent+child inputs into a dedicated decomposition family rather than generic structural reporting.",
            "The current implementation preserves detected magnetic metadata for downstream decomposition tooling.",
        ],
        "incommensurate_magnetic_superspace_analysis": [
            "Routes incommensurate magnetic cases into a dedicated superspace family.",
            "This wrapper records superspace evidence explicitly so a conventional 3D structural workflow is not selected accidentally.",
        ],
        "mixed_magnetic_structural_coupled_analysis": [
            "Routes mixed magnetic and structural distortions into a coupled-analysis family.",
            "This wrapper keeps both magnetic and structural metadata visible for future joint OP analysis.",
        ],
    }[workflow]

    return {
        "label": label,
        "workflow": workflow,
        "title": WORKFLOW_TITLES[workflow],
        "inputs": {
            "parent_cif": str(parent_path) if parent_path else None,
            "child_cif": str(child_path) if child_path else None,
            "distortion_file": str(distortion_path) if distortion_path else None,
        },
        "metadata": metadata,
        "notes": notes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a dedicated magnetic/superspace workflow scaffold.")
    parser.add_argument("--workflow", required=True, choices=sorted(WORKFLOW_TITLES))
    parser.add_argument("--label", default="magnetic_run")
    parser.add_argument("--parent-cif")
    parser.add_argument("--child-cif")
    parser.add_argument("--distortion-file")
    args = parser.parse_args()

    run_dir = MAGNETIC_RUNS_DIR / args.label
    run_dir.mkdir(parents=True, exist_ok=True)
    summary = build_summary(
        workflow=args.workflow,
        label=args.label,
        parent_cif=args.parent_cif,
        child_cif=args.child_cif,
        distortion_file=args.distortion_file,
    )
    out = run_dir / "magnetic_workflow_report.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

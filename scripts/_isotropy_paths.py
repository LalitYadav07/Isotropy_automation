from __future__ import annotations

import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
VENDOR_DIR = PROJECT_ROOT / "vendor" / "isobyu"
INPUTS_DIR = PROJECT_ROOT / "inputs"
TUTORIALS_DIR = PROJECT_ROOT / "isotutorials"
OUTPUTS_DIR = PROJECT_ROOT / "outputs"
SUBGROUPS_DIR = OUTPUTS_DIR / "subgroups"
PLOTS_DIR = OUTPUTS_DIR / "diffraction_plots"
REPORT_PATH = OUTPUTS_DIR / "scientific_report.md"
TUTORIAL_REPORTS_DIR = OUTPUTS_DIR / "tutorial_reports"
DISCOVERY_RUNS_DIR = OUTPUTS_DIR / "discovery_runs"
SCRATCH_DIR = PROJECT_ROOT / "scratch"


def configure_environment() -> None:
    for path in (SUBGROUPS_DIR, PLOTS_DIR, TUTORIAL_REPORTS_DIR, DISCOVERY_RUNS_DIR, SCRATCH_DIR):
        path.mkdir(parents=True, exist_ok=True)

    os.environ["ISODATA"] = f"{VENDOR_DIR}{os.sep}"
    current_path = os.environ.get("PATH", "")
    vendor_str = str(VENDOR_DIR)
    if current_path:
        path_parts = current_path.split(os.pathsep)
        if vendor_str not in path_parts:
            os.environ["PATH"] = os.pathsep.join([vendor_str, current_path])
    else:
        os.environ["PATH"] = vendor_str


def resolve_input_path(path_str: str) -> Path:
    path = Path(path_str)
    if path.is_absolute():
        return path

    cwd_path = Path.cwd() / path
    if cwd_path.exists():
        return cwd_path.resolve()

    return (INPUTS_DIR / path).resolve()

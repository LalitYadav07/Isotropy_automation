from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from _isotropy_paths import SCRATCH_DIR, configure_environment


configure_environment()


@dataclass
class StandardizationResult:
    standardized_path: Path
    structure: Structure
    method: str
    raw_space_group_number: int
    raw_space_group_symbol: str
    standardized_space_group_number: int
    standardized_space_group_symbol: str


def _safe_slug(path: Path) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", path.stem)


def _run_findsym_from_cif(source_cif: Path, workdir: Path, tolerance: float) -> Path | None:
    temp_input = workdir / f"{_safe_slug(source_cif)}_findsym.in"
    generated_cif = workdir / "findsym.cif"
    generated_log = workdir / "findsym.log"

    for artifact in (temp_input, generated_cif, generated_log):
        if artifact.exists():
            artifact.unlink()

    with temp_input.open("w", encoding="utf-8") as handle:
        subprocess.run(
            ["findsym_cifinput", str(source_cif)],
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
            cwd=workdir,
            check=False,
        )

    subprocess.run(
        ["findsym", str(temp_input)],
        input=f"\n{tolerance}\n",
        capture_output=True,
        text=True,
        cwd=workdir,
        check=False,
    )

    if generated_cif.exists():
        return generated_cif
    return None


def standardize_cif_with_fallback(
    source_cif: Path,
    output_cif: Path,
    *,
    symprec: float = 1e-3,
    findsym_tolerance: float = 1e-3,
) -> StandardizationResult:
    raw_structure = CifParser(str(source_cif)).parse_structures(primitive=False)[0]
    raw_analyzer = SpacegroupAnalyzer(raw_structure, symprec=symprec)

    findsym_workspace = SCRATCH_DIR / "findsym_standardize"
    findsym_workspace.mkdir(parents=True, exist_ok=True)
    generated_cif = _run_findsym_from_cif(source_cif, findsym_workspace, findsym_tolerance)

    standardized_structure: Structure
    method: str
    if generated_cif is not None:
        standardized_structure = CifParser(str(generated_cif)).parse_structures(primitive=False)[0]
        method = "findsym"
        shutil.copyfile(generated_cif, output_cif)
    else:
        standardized_structure = raw_analyzer.get_conventional_standard_structure()
        CifWriter(standardized_structure).write_file(output_cif)
        method = "pymatgen_conventional_standard"

    analyzer = SpacegroupAnalyzer(standardized_structure, symprec=symprec)
    return StandardizationResult(
        standardized_path=output_cif,
        structure=standardized_structure,
        method=method,
        raw_space_group_number=raw_analyzer.get_space_group_number(),
        raw_space_group_symbol=raw_analyzer.get_space_group_symbol(),
        standardized_space_group_number=analyzer.get_space_group_number(),
        standardized_space_group_symbol=analyzer.get_space_group_symbol(),
    )

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from _isotropy_paths import OUTPUTS_DIR, resolve_input_path
from findsym_tools import standardize_cif_with_fallback


RECONSTRUCTIVE_DIR = OUTPUTS_DIR / "reconstructive_runs"
RECONSTRUCTIVE_DIR.mkdir(parents=True, exist_ok=True)


def _wyckoff_rows(structure) -> list[tuple[str, str, list[float]]]:
    analyzer = SpacegroupAnalyzer(structure, symprec=1e-2)
    symm = analyzer.get_symmetrized_structure()
    rows = []
    for site, wyckoff in zip(symm, symm.wyckoff_symbols, strict=False):
        rows.append((site.species_string, wyckoff[-1].lower(), [float(value) for value in site.frac_coords]))
    unique_rows: list[tuple[str, str, list[float]]] = []
    seen = set()
    for species, wyckoff, coords in rows:
        key = (species, wyckoff, tuple(round(value, 6) for value in coords))
        if key in seen:
            continue
        seen.add(key)
        unique_rows.append((species, wyckoff, coords))
    return unique_rows


def _lattice_line(structure) -> str:
    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles
    return f"{a:.6f} {b:.6f} {c:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f}"


def build_comsubs_input(parent_cif: Path, child_cif: Path, label: str) -> tuple[str, dict[str, object]]:
    run_dir = RECONSTRUCTIVE_DIR / label
    run_dir.mkdir(parents=True, exist_ok=True)

    parent_std = standardize_cif_with_fallback(parent_cif, run_dir / "parent_standardized.cif")
    child_std = standardize_cif_with_fallback(child_cif, run_dir / "child_standardized.cif")

    parent_rows = _wyckoff_rows(parent_std.structure)
    child_rows = _wyckoff_rows(child_std.structure)

    lines = [
        f"Reconstructive comparison for {parent_cif.name} -> {child_cif.name}",
        str(parent_std.standardized_space_group_number),
        _lattice_line(parent_std.structure),
        str(len(parent_rows)),
    ]
    lines.extend(
        f"{species} {wyckoff} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}"
        for species, wyckoff, coords in parent_rows
    )
    lines.extend(
        [
            str(child_std.standardized_space_group_number),
            _lattice_line(child_std.structure),
            str(len(child_rows)),
        ]
    )
    lines.extend(
        f"{species} {wyckoff} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}"
        for species, wyckoff, coords in child_rows
    )
    lines.extend(
        [
            "size: 1 16",
            "strain: 0.5 2.0",
            "neighbor: 0.0",
            "shuffle: 5.0",
        ]
    )

    metadata = {
        "parent_standardization_method": parent_std.method,
        "child_standardization_method": child_std.method,
        "parent_space_group": f"{parent_std.standardized_space_group_number} {parent_std.standardized_space_group_symbol}",
        "child_space_group": f"{child_std.standardized_space_group_number} {child_std.standardized_space_group_symbol}",
        "parent_wyckoff_count": len(parent_rows),
        "child_wyckoff_count": len(child_rows),
    }
    return "\n".join(lines) + "\n", metadata


def parse_comsubs_candidates(text: str) -> list[dict[str, object]]:
    candidates = []
    pattern = re.compile(r"^\s*(\d+)\s+([A-Za-z0-9/_-]+)\s+size\s*=\s*(\d+)", re.IGNORECASE)
    for line in text.splitlines():
        match = pattern.search(line)
        if not match:
            continue
        candidates.append(
            {
                "space_group_number": int(match.group(1)),
                "space_group_symbol": match.group(2),
                "cell_size": int(match.group(3)),
                "raw_line": line.strip(),
            }
        )
    return candidates


def run_comsubs(parent_cif: str, child_cif: str, label: str) -> dict[str, object]:
    parent_path = resolve_input_path(parent_cif)
    child_path = resolve_input_path(child_cif)
    run_dir = RECONSTRUCTIVE_DIR / label
    run_dir.mkdir(parents=True, exist_ok=True)

    comsubs_input, metadata = build_comsubs_input(parent_path, child_path, label)
    input_path = run_dir / "comsubs.in"
    input_path.write_text(comsubs_input, encoding="utf-8")

    result = subprocess.run(
        ["comsubs"],
        input=comsubs_input,
        text=True,
        capture_output=True,
        cwd=run_dir,
        check=False,
    )
    stdout_path = run_dir / "comsubs.out"
    stdout_path.write_text(result.stdout, encoding="utf-8")
    if result.stderr:
        (run_dir / "comsubs.err").write_text(result.stderr, encoding="utf-8")

    candidates = parse_comsubs_candidates(result.stdout)
    report = {
        "label": label,
        "parent_cif": str(parent_path),
        "child_cif": str(child_path),
        "metadata": metadata,
        "return_code": result.returncode,
        "candidate_count": len(candidates),
        "candidates": candidates,
        "artifacts": {
            "input": str(input_path),
            "stdout": str(stdout_path),
        },
        "notes": [
            "This wrapper treats comsubs as a reconstructive/common-subgroup workflow, not as an extension of parent-only irrep discovery.",
            "Candidate parsing is intentionally conservative; the raw comsubs output is preserved alongside the machine summary.",
        ],
    }
    (run_dir / "reconstructive_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    markdown_lines = [
        "# Reconstructive / Common-Subgroup Report",
        "",
        f"- Parent CIF: `{report['parent_cif']}`",
        f"- Child CIF: `{report['child_cif']}`",
        f"- Parent SG: `{metadata['parent_space_group']}`",
        f"- Child SG: `{metadata['child_space_group']}`",
        f"- Parent standardization: `{metadata['parent_standardization_method']}`",
        f"- Child standardization: `{metadata['child_standardization_method']}`",
        f"- Candidate count: `{report['candidate_count']}`",
        f"- `comsubs` return code: `{report['return_code']}`",
        "",
        "## Notes",
    ]
    markdown_lines.extend(f"- {note}" for note in report["notes"])
    markdown_lines.extend(
        [
            "",
            "## Candidates",
        ]
    )
    if candidates:
        markdown_lines.extend(
            f"- `{candidate['space_group_number']} {candidate['space_group_symbol']}` size `{candidate['cell_size']}` from `{candidate['raw_line']}`"
            for candidate in candidates
        )
    else:
        markdown_lines.append("- No conservative machine-parsed candidates were extracted; inspect the raw `comsubs.out` artifact directly.")
    (run_dir / "reconstructive_report.md").write_text("\n".join(markdown_lines) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Search reconstructive/common-subgroup pathways with comsubs.")
    parser.add_argument("parent_cif")
    parser.add_argument("child_cif")
    parser.add_argument("--label")
    args = parser.parse_args()

    label = args.label or f"{Path(args.parent_cif).stem}_to_{Path(args.child_cif).stem}"
    report = run_comsubs(args.parent_cif, args.child_cif, label)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

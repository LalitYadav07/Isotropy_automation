from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

try:
    from _isotropy_paths import OUTPUTS_DIR, configure_environment, resolve_input_path
    from findsym_tools import standardize_cif_with_fallback
    from isodistort_distortion_parser import DistortionReportData, parse_distortion_file
except ModuleNotFoundError:
    from isotropy_project.scripts._isotropy_paths import OUTPUTS_DIR, configure_environment, resolve_input_path
    from isotropy_project.scripts.findsym_tools import standardize_cif_with_fallback
    from isotropy_project.scripts.isodistort_distortion_parser import DistortionReportData, parse_distortion_file


configure_environment()

MAGNETIC_CHILD_RUNS_DIR = OUTPUTS_DIR / "magnetic_child_runs"
MAGNETIC_CHILD_RUNS_DIR.mkdir(parents=True, exist_ok=True)


def parse_order_parameter_details(order_parameter: str) -> dict[str, str]:
    details: dict[str, str] = {"raw": order_parameter}
    lead = order_parameter.split(", basis=", 1)[0].strip()
    basis_match = re.search(r"basis=\{([^}]*)\}", order_parameter)
    origin_match = re.search(r"origin=\(([^)]*)\)", order_parameter)
    size_match = re.search(r"\bs=([^,]+)", order_parameter)
    index_match = re.search(r"\bi=([^,]+)", order_parameter)

    if lead.startswith("Subgroup: "):
        details["subgroup_selection"] = lead.removeprefix("Subgroup: ").strip()
    else:
        parts = lead.split(" ", 2)
        if parts:
            details["opd_label"] = parts[0]
        if len(parts) > 1:
            details["opd_vector"] = parts[1]
        if len(parts) > 2:
            details["subgroup_selection"] = parts[2]

    if basis_match:
        details["basis"] = basis_match.group(1)
    if origin_match:
        details["origin"] = origin_match.group(1)
    if size_match:
        details["size_index"] = size_match.group(1).strip()
    if index_match:
        details["symmetry_index"] = index_match.group(1).strip()
    return details


def _structure_and_sg(path: Path) -> tuple[object, str]:
    structure = CifParser(str(path)).parse_structures(primitive=False)[0]
    analyzer = SpacegroupAnalyzer(structure, symprec=1e-2)
    return structure, f"{analyzer.get_space_group_number()} {analyzer.get_space_group_symbol()}"


def _composition_ratio(parent_structure, child_structure) -> float | None:
    parent_atoms = parent_structure.composition.num_atoms
    child_atoms = child_structure.composition.num_atoms
    if parent_atoms <= 0:
        return None
    return child_atoms / parent_atoms


def _channel_family_label(family: str) -> str:
    return {
        "magnetic": "magnetic primary-order candidate",
        "displacive": "secondary displacive",
        "ordering": "secondary occupational/order-disorder",
        "rotational": "secondary rotational",
        "strain": "secondary strain",
        "ellipsoidal": "secondary ellipsoidal",
    }.get(family, family)


def _collect_channel_rows(data: DistortionReportData) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    magnetic_rss = float(data.mode_summaries.get("magnetic", {}).get("rss", 0.0) or 0.0)
    for family in ("magnetic", "displacive", "ordering", "rotational", "strain", "ellipsoidal"):
        detail = data.mode_summaries.get(family)
        if not detail:
            continue
        rss = float(detail.get("rss", 0.0) or 0.0)
        rows.append(
            {
                "family": family,
                "label": _channel_family_label(family),
                "rss": rss,
                "nonzero_mode_count": int(detail.get("nonzero_mode_count", 0) or 0),
                "max_abs": float(detail.get("max_abs", 0.0) or 0.0),
                "relative_to_magnetic_rss": (rss / magnetic_rss) if magnetic_rss > 0 else None,
            }
        )
    return rows


def _identify_primary_magnetic_contributors(data: DistortionReportData) -> list[dict[str, object]]:
    magnetic = data.mode_summaries.get("magnetic", {})
    per_irrep = magnetic.get("per_irrep", {}) if isinstance(magnetic, dict) else {}
    candidates: list[dict[str, object]] = []

    primary_by_index = {item.get("index"): item for item in data.primary_irreps}
    for irrep_index_str, summary in per_irrep.items():
        rss = float(summary.get("rss", 0.0) or 0.0)
        if rss <= 1e-9:
            continue
        irrep_index = int(irrep_index_str)
        primary_entry = primary_by_index.get(irrep_index)
        candidates.append(
            {
                "irrep_slot": irrep_index,
                "irrep_label": primary_entry.get("irrep") if primary_entry else None,
                "k_vector": primary_entry.get("k_vector") if primary_entry else None,
                "rss": rss,
                "nonzero_mode_count": int(summary.get("nonzero_count", 0.0) or 0),
                "max_abs": float(summary.get("max_abs", 0.0) or 0.0),
                "export_label_status": "explicit" if primary_entry and primary_entry.get("irrep") else "slot_only",
            }
        )
    candidates.sort(key=lambda item: item["rss"], reverse=True)
    total = sum(item["rss"] ** 2 for item in candidates)
    for item in candidates:
        item["fraction_of_magnetic_rss"] = (item["rss"] ** 2 / total) if total > 0 else None
    return candidates


def build_report(parent_cif: str, child_cif: str, distortion_file: str, label: str) -> dict[str, object]:
    parent_path = resolve_input_path(parent_cif)
    child_path = resolve_input_path(child_cif)
    distortion_path = resolve_input_path(distortion_file)

    run_dir = MAGNETIC_CHILD_RUNS_DIR / label
    run_dir.mkdir(parents=True, exist_ok=True)

    parent_std = standardize_cif_with_fallback(parent_path, run_dir / "parent_standardized.cif")
    child_std = standardize_cif_with_fallback(child_path, run_dir / "child_standardized.cif")
    parent_structure, parent_sg = _structure_and_sg(parent_path)
    child_structure, child_sg = _structure_and_sg(child_path)
    matcher = StructureMatcher(primitive_cell=False, scale=True, attempt_supercell=True)
    same_framework = bool(matcher.fit(parent_std.structure, child_std.structure))

    distortion = parse_distortion_file(distortion_path)
    opd = parse_order_parameter_details(distortion.order_parameter)
    channels = _collect_channel_rows(distortion)
    primary_magnetic = _identify_primary_magnetic_contributors(distortion)
    secondary_channels = [row for row in channels if row["family"] != "magnetic" and row["rss"] > 1e-9]
    secondary_channels.sort(key=lambda item: item["rss"], reverse=True)
    dominant_secondary = secondary_channels[0] if secondary_channels else None

    composition_ratio = _composition_ratio(parent_std.structure, child_std.structure)
    magnetic_rss = float(distortion.mode_summaries.get("magnetic", {}).get("rss", 0.0) or 0.0)
    secondary_rss = math.sqrt(sum(row["rss"] ** 2 for row in secondary_channels)) if secondary_channels else 0.0

    report = {
        "label": label,
        "inputs": {
            "parent_cif": str(parent_path),
            "child_cif": str(child_path),
            "distortion_file": str(distortion_path),
        },
        "setting_alignment": {
            "parent_raw_space_group": parent_sg,
            "child_raw_space_group": child_sg,
            "parent_standardized_space_group": f"{parent_std.standardized_space_group_number} {parent_std.standardized_space_group_symbol}",
            "child_standardized_space_group": f"{child_std.standardized_space_group_number} {child_std.standardized_space_group_symbol}",
            "parent_standardization_method": parent_std.method,
            "child_standardization_method": child_std.method,
            "structure_matcher_fit_after_standardization": same_framework,
            "child_to_parent_atom_count_ratio": composition_ratio,
            "parent_formula": parent_std.structure.composition.reduced_formula,
            "child_formula": child_std.structure.composition.reduced_formula,
        },
        "opd_summary": {
            "order_parameter": distortion.order_parameter,
            "parsed": opd,
            "subgroup_string": distortion.subgroup_string,
        },
        "primary_magnetic_order": {
            "contributors": primary_magnetic,
            "dominant": primary_magnetic[0] if primary_magnetic else None,
            "note": (
                "Exact irrep labels are available only when the distortion export stores them explicitly. "
                "Mode-decomposition exports often preserve only irrep slots, so the report falls back to those slots when needed."
            ),
        },
        "channel_amplitudes": channels,
        "secondary_channels": secondary_channels,
        "physical_summary": {
            "magnetic_rss_amplitude": magnetic_rss,
            "secondary_combined_rss_amplitude": secondary_rss,
            "magnetic_to_secondary_ratio": (magnetic_rss / secondary_rss) if secondary_rss > 0 else None,
            "dominant_secondary_channel": dominant_secondary,
            "interpretation": [
                "Treat magnetic RSS amplitude as the size of the populated magnetic symmetry-mode manifold in the exported decomposition, not as an ordered moment on a single site.",
                "Treat non-magnetic RSS amplitudes as secondary responses living in the same subgroup setting: displacive, occupational, rotational, and strain channels can all be symmetry-allowed once the magnetic primary order is selected.",
                "If the dominant magnetic contributor is the only appreciably populated magnetic irrep slot, that slot is the best automated proxy for the primary magnetic order parameter in this child structure.",
            ],
        },
        "distortion_checks": distortion.internal_checks,
        "classification": distortion.classification,
    }

    (run_dir / "magnetic_child_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = [
        "# Magnetic Parent/Child Decomposition Report",
        "",
        "## Inputs",
        f"- Parent CIF: `{parent_path}`",
        f"- Magnetic child CIF: `{child_path}`",
        f"- Distortion file: `{distortion_path}`",
        "",
        "## Setting alignment",
        f"- Parent raw space group: `{parent_sg}`",
        f"- Child raw space group: `{child_sg}`",
        f"- Parent standardized setting: `{report['setting_alignment']['parent_standardized_space_group']}` via `{parent_std.method}`",
        f"- Child standardized setting: `{report['setting_alignment']['child_standardized_space_group']}` via `{child_std.method}`",
        f"- StructureMatcher fit after standardization: `{same_framework}`",
    ]
    if composition_ratio is not None:
        lines.append(f"- Child/parent atom-count ratio: `{composition_ratio:.3f}`")
    lines.extend(
        [
            "",
            "## Primary magnetic order",
        ]
    )
    if primary_magnetic:
        for item in primary_magnetic:
            label_text = item['irrep_label'] or f"slot {item['irrep_slot']}"
            fraction = item.get('fraction_of_magnetic_rss')
            fraction_text = f", fraction={fraction:.3f}" if fraction is not None else ""
            lines.append(
                f"- `{label_text}`: rss `{item['rss']:.6f}`, nonzero modes `{item['nonzero_mode_count']}`, max |coef| `{item['max_abs']:.6f}`{fraction_text}"
            )
    else:
        lines.append("- No nonzero magnetic mode amplitudes were detected in the supplied distortion file.")
    lines.extend(
        [
            "",
            "## OPD and subgroup",
            f"- Stored order parameter: `{distortion.order_parameter}`",
            f"- Parsed subgroup selection: `{opd.get('subgroup_selection', 'not explicit')}`",
            f"- Parsed basis: `{opd.get('basis', 'not explicit')}`",
            f"- Parsed origin: `{opd.get('origin', 'not explicit')}`",
            "",
            "## Secondary channels",
        ]
    )
    if secondary_channels:
        for row in secondary_channels:
            rel = row['relative_to_magnetic_rss']
            rel_text = f", relative to magnetic rss `{rel:.3f}`" if rel is not None else ""
            lines.append(
                f"- `{row['label']}`: rss `{row['rss']:.6f}`, nonzero modes `{row['nonzero_mode_count']}`, max |coef| `{row['max_abs']:.6f}`{rel_text}"
            )
    else:
        lines.append("- No non-magnetic secondary channels carry nonzero amplitude in this export.")
    lines.extend(
        [
            "",
            "## Physical interpretation",
            f"- Magnetic RSS amplitude: `{magnetic_rss:.6f}`",
            f"- Combined secondary RSS amplitude: `{secondary_rss:.6f}`",
            f"- Magnetic/secondary ratio: `{report['physical_summary']['magnetic_to_secondary_ratio']}`",
        ]
    )
    if dominant_secondary:
        lines.append(
            f"- Dominant secondary channel: `{dominant_secondary['label']}` with rss `{dominant_secondary['rss']:.6f}`"
        )
    lines.extend(["", "## Notes"])
    lines.extend(f"- {note}" for note in report["physical_summary"]["interpretation"])
    (run_dir / "magnetic_child_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze a known magnetic child structure relative to a known parent using an exported ISODISTORT distortion file."
    )
    parser.add_argument("parent_cif")
    parser.add_argument("child_cif")
    parser.add_argument("distortion_file")
    parser.add_argument("--label")
    args = parser.parse_args()

    label = args.label or f"{Path(args.parent_cif).stem}_to_{Path(args.child_cif).stem}_magnetic"
    report = build_report(args.parent_cif, args.child_cif, args.distortion_file, label)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

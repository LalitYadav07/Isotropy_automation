from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
import math


@dataclass
class DistortionReportData:
    source_path: str
    slug: str
    title: str
    parent_string: str
    parent_formula_estimate: str
    parent_setting: str | None
    lattice_string: str
    order_parameter: str
    subgroup_string: str | None
    wyckoff_strings: list[str]
    k_vectors: list[str]
    irreps: list[str]
    distortion_irrep_numbers: list[int]
    include_flags: dict[str, list[str] | str]
    mode_counts: dict[str, int]
    coefficient_summaries: dict[str, dict[str, float]]
    mode_summaries: dict[str, dict[str, object]]
    primary_irreps: list[dict[str, object]]
    modes_irrep_numbers: list[int]
    irrep_count: int
    modulation_count: int
    subgroup_cell_size: int
    classification: list[str]
    associated_files: list[str]
    tags: dict[str, str]
    child_basic_structure: dict[str, object] | None
    internal_checks: dict[str, object]


def _parse_sections(text: str) -> dict[str, dict[str, str]]:
    sections: dict[str, dict[str, str]] = {}
    current_section: str | None = None
    current_tag: str | None = None
    buffer: list[str] = []

    def flush() -> None:
        nonlocal buffer
        if current_section and buffer:
            key = current_tag or "__body__"
            content = "\n".join(buffer).strip()
            if content:
                section = sections.setdefault(current_section, {})
                existing = section.get(key)
                section[key] = f"{existing}\n{content}".strip() if existing else content
        buffer = []

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith("!begin "):
            flush()
            current_section = line.split(maxsplit=1)[1]
            current_tag = None
            continue
        if line.startswith("!end "):
            flush()
            current_section = None
            current_tag = None
            continue
        if current_section is None:
            continue
        if line.startswith("!"):
            flush()
            current_tag = line[1:].strip()
            continue
        if line.startswith("#"):
            continue
        buffer.append(line)
    flush()
    return sections


def _lines(value: str | None) -> list[str]:
    if not value:
        return []
    return [line.strip() for line in value.splitlines() if line.strip()]


def _numbers(value: str | None) -> list[float]:
    if not value:
        return []
    values: list[float] = []
    for token in value.replace("\n", " ").split():
        try:
            values.append(float(token))
        except ValueError:
            continue
    return values


def _ints(value: str | None) -> list[int]:
    return [int(round(number)) for number in _numbers(value)]


def _count_nonzero(values: list[float], tol: float = 1e-9) -> int:
    return sum(1 for value in values if abs(value) > tol)


def _rss(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def _bool_tokens(value: str | None) -> list[str]:
    if not value:
        return []
    return [token for token in value.replace("\n", " ").split() if token in {"T", "F"}]


def _classification(
    include_flags: dict[str, list[str] | str],
    modulation_count: int,
    distortion_irrep_count: int,
) -> list[str]:
    classes: list[str] = []
    if "T" in include_flags.get("includeMagnetic", []):
        classes.append("magnetic")
    if "T" in include_flags.get("includeDisplacive", []):
        classes.append("displacive")
    if "T" in include_flags.get("includeOrdering", []):
        classes.append("occupational")
    if "T" in include_flags.get("includeRotational", []):
        classes.append("rotational")
    if include_flags.get("includeStrain") == "T":
        classes.append("strain")
    classes.append("incommensurate/modulated" if modulation_count > 0 else "commensurate")
    classes.append("mode decomposition" if distortion_irrep_count == 0 else "generated distortion")
    return classes


def _associated_files(path: Path) -> list[str]:
    prefix = path.stem.replace("-distortion", "")
    results: list[str] = []
    for candidate in sorted(path.parent.iterdir()):
        if candidate.name == path.name:
            continue
        if candidate.stem.startswith(prefix):
            results.append(candidate.name)
    return results


def _format_formula(parts: list[tuple[str, float]]) -> str:
    tokens: list[str] = []
    for symbol, count in parts:
        rounded = round(count)
        if abs(count - rounded) < 1e-8:
            count_text = "" if rounded == 1 else str(int(rounded))
        else:
            count_text = f"{count:.3f}".rstrip("0").rstrip(".")
        tokens.append(f"{symbol}{count_text}")
    return "".join(tokens)


def _reduced_integer_counts(parts: list[tuple[str, float]]) -> list[tuple[str, float]]:
    rounded = [int(round(count)) for _, count in parts]
    if not all(abs(count - rounded_value) < 1e-8 for (_, count), rounded_value in zip(parts, rounded, strict=False)):
        return parts

    gcd_value = 0
    for value in rounded:
        gcd_value = math.gcd(gcd_value, value)
    if gcd_value <= 1:
        return parts
    return [(symbol, count / gcd_value) for symbol, count in parts]


def _parent_formula_estimate(distortion: dict[str, str]) -> str:
    symbols = _lines(distortion.get("wyckoffAtom"))
    multiplicities = []
    for wyckoff in _lines(distortion.get("wyckoffString")):
        match = wyckoff.split(maxsplit=1)[0]
        digits = "".join(character for character in match if character.isdigit())
        multiplicities.append(int(digits) if digits else 1)
    occupations = _numbers(distortion.get("wyckoffOccupation"))
    if not symbols or not multiplicities:
        return ""
    if len(occupations) < len(symbols):
        occupations.extend([1.0] * (len(symbols) - len(occupations)))

    totals: dict[str, float] = defaultdict(float)
    order: list[str] = []
    for index, symbol in enumerate(symbols):
        multiplicity = multiplicities[index] if index < len(multiplicities) else 1
        occupation = occupations[index] if index < len(occupations) else 1.0
        if symbol not in totals:
            order.append(symbol)
        totals[symbol] += float(multiplicity) * float(occupation)
    return _format_formula(_reduced_integer_counts([(symbol, totals[symbol]) for symbol in order]))


def _trim_site_coefficients(values: list[float], counts: list[int], max_per_site: int) -> list[float]:
    if not counts or max_per_site <= 0:
        return values

    trimmed: list[float] = []
    cursor = 0
    for count in counts:
        trimmed.extend(values[cursor : cursor + count])
        cursor += max_per_site
    return trimmed


def _summarize_coefficients(values: list[float]) -> dict[str, float]:
    return {
        "count": float(len(values)),
        "nonzero_count": float(_count_nonzero(values)),
        "rss": _rss(values),
        "max_abs": max((abs(value) for value in values), default=0.0),
    }


def _mode_family_summary(
    family: str,
    distortion: dict[str, str],
    modes: dict[str, str],
) -> tuple[dict[str, float], dict[str, object]]:
    cap = family[0].upper() + family[1:]
    count_tag = f"{family}ModesCount"
    coeff_tag = f"{family}ModesCoef"
    irrep_tag = f"{family}ModeIrrep"

    if family == "strain":
        counts = _ints(modes.get(count_tag))
        total_mode_count = int(sum(counts))
        coefficients = _numbers(distortion.get(coeff_tag))[:total_mode_count]
        irrep_assignments = _ints(modes.get(irrep_tag))
    else:
        counts = _ints(modes.get(count_tag))
        total_mode_count = int(sum(counts))
        max_per_site = int(sum(_numbers(distortion.get(f"max{cap}Modes")) or _numbers(modes.get(f"max{cap}Modes")) or [0]))
        raw_coefficients = _numbers(distortion.get(coeff_tag))
        coefficients = _trim_site_coefficients(raw_coefficients, counts, max_per_site)
        irrep_assignments = _ints(modes.get(irrep_tag))

    per_irrep: dict[str, dict[str, float]] = {}
    for coefficient, irrep_index in zip(coefficients, irrep_assignments, strict=False):
        bucket = per_irrep.setdefault(
            str(irrep_index),
            {"count": 0.0, "nonzero_count": 0.0, "rss": 0.0, "max_abs": 0.0},
        )
        bucket["count"] += 1.0
        if abs(coefficient) > 1e-9:
            bucket["nonzero_count"] += 1.0
            bucket["rss"] += coefficient * coefficient
            bucket["max_abs"] = max(bucket["max_abs"], abs(coefficient))

    for bucket in per_irrep.values():
        bucket["rss"] = math.sqrt(bucket["rss"])

    summary = _summarize_coefficients(coefficients)
    summary["total_mode_count"] = float(total_mode_count)
    summary["assignment_count"] = float(len(irrep_assignments))

    detail = {
        "total_mode_count": total_mode_count,
        "coefficient_count": len(coefficients),
        "assignment_count": len(irrep_assignments),
        "nonzero_mode_count": _count_nonzero(coefficients),
        "rss": summary["rss"],
        "max_abs": summary["max_abs"],
        "consistency_ok": total_mode_count == len(coefficients) == len(irrep_assignments),
        "per_irrep": per_irrep,
    }
    return summary, detail


def _parse_atoms_body(body: str | None) -> dict[str, object] | None:
    lines = _lines(body)
    if not lines:
        return None

    lattice_parameters = _numbers(lines[0])[:6]
    atoms: list[dict[str, object]] = []
    index = 1
    while index + 1 < len(lines):
        label_tokens = lines[index].split()
        if len(label_tokens) >= 2:
            label, species = label_tokens[0], label_tokens[1]
        elif label_tokens:
            label = species = label_tokens[0]
        else:
            index += 2
            continue
        coords = _numbers(lines[index + 1])[:3]
        atoms.append(
            {
                "label": label,
                "species": species,
                "frac_coords": coords,
            }
        )
        index += 2

    return {
        "lattice_parameters": lattice_parameters,
        "atom_count": len(atoms),
        "atoms": atoms,
    }


def _primary_irrep_summaries(
    distortion_irrep_numbers: list[int],
    irreps: list[str],
    k_vectors: list[str],
    mode_summaries: dict[str, dict[str, object]],
    modes_irrep_numbers: list[int],
) -> tuple[list[dict[str, object]], bool]:
    if not distortion_irrep_numbers:
        return [], True

    alignment_ok = modes_irrep_numbers[: len(distortion_irrep_numbers)] == distortion_irrep_numbers
    results: list[dict[str, object]] = []
    for index, (irrep_label, k_vector, irrep_number) in enumerate(
        zip_longest(irreps, k_vectors, distortion_irrep_numbers, fillvalue=""),
        start=1,
    ):
        category_summaries: dict[str, dict[str, float]] = {}
        if alignment_ok:
            for family, summary in mode_summaries.items():
                category_summaries[family] = summary["per_irrep"].get(
                    str(index),
                    {"count": 0.0, "nonzero_count": 0.0, "rss": 0.0, "max_abs": 0.0},
                )
        results.append(
            {
                "index": index,
                "irrep": irrep_label,
                "k_vector": k_vector,
                "irrep_number": irrep_number,
                "category_summaries": category_summaries,
            }
        )
    return results, alignment_ok


def parse_distortion_file(path: Path) -> DistortionReportData:
    text = path.read_text(encoding="utf-8")
    sections = _parse_sections(text)
    distortion = sections.get("distortionFile", {})
    modes = sections.get("modesFile", {})
    atoms = sections.get("atomsFile", {})

    include_flags = {
        "includeDisplacive": _bool_tokens(distortion.get("includeDisplacive")),
        "includeOrdering": _bool_tokens(distortion.get("includeOrdering")),
        "includeMagnetic": _bool_tokens(distortion.get("includeMagnetic")),
        "includeRotational": _bool_tokens(distortion.get("includeRotational")),
        "includeEllipsoidal": _bool_tokens(distortion.get("includeEllipsoidal")),
        "includeStrain": (distortion.get("includeStrain") or "").strip(),
    }

    mode_counts = {
        "displacive": int(sum(_numbers(modes.get("displaciveModesCount")))),
        "magnetic": int(sum(_numbers(modes.get("magneticModesCount")))),
        "rotational": int(sum(_numbers(modes.get("rotationalModesCount")))),
        "ellipsoidal": int(sum(_numbers(modes.get("ellipsoidalModesCount")))),
        "ordering": int(sum(_numbers(modes.get("orderingModesCount")))),
        "strain": int(sum(_numbers(modes.get("strainModesCount")))),
    }

    coefficient_summaries: dict[str, dict[str, float]] = {}
    mode_summaries: dict[str, dict[str, object]] = {}
    for family in ("displacive", "magnetic", "ordering", "rotational", "ellipsoidal", "strain"):
        summary, detail = _mode_family_summary(family, distortion, modes)
        coefficient_summaries[f"{family}ModesCoef"] = summary
        mode_summaries[family] = detail

    modulation_count = int(sum(_numbers(modes.get("modCount")) or _numbers(distortion.get("modCount")) or [0]))
    subgroup_cell_size = int(
        sum(_numbers(modes.get("subgroupCellSize")) or _numbers(distortion.get("subgroupCellSize")) or [0])
    )
    modes_irrep_numbers = _ints(modes.get("irrepNumber"))
    distortion_irrep_numbers = _ints(distortion.get("irrepNumber"))
    irreps = _lines(distortion.get("irrepString"))
    k_vectors = _lines(distortion.get("kvecString"))
    primary_irreps, primary_irrep_alignment = _primary_irrep_summaries(
        distortion_irrep_numbers,
        irreps,
        k_vectors,
        mode_summaries,
        modes_irrep_numbers,
    )

    slug = path.stem[:-11] if path.stem.endswith("-distortion") else path.stem
    title = slug.replace("-", " ")
    child_basic_structure = _parse_atoms_body(atoms.get("__body__"))
    associated_files = _associated_files(path)
    missing_associated_files = [
        filename for filename in associated_files if not (path.parent / filename).exists()
    ]
    coefficient_consistency_ok = all(summary["consistency_ok"] for summary in mode_summaries.values())

    return DistortionReportData(
        source_path=str(path),
        slug=slug,
        title=title,
        parent_string=(distortion.get("parentString") or "").strip(),
        parent_formula_estimate=_parent_formula_estimate(distortion),
        parent_setting=(distortion.get("parentSettingString") or "").strip() or None,
        lattice_string=(distortion.get("lattParamString") or "").strip(),
        order_parameter=(distortion.get("orderParamString") or "").strip(),
        subgroup_string=(distortion.get("subgroupString") or "").strip() or None,
        wyckoff_strings=_lines(distortion.get("wyckoffString")),
        k_vectors=k_vectors,
        irreps=irreps,
        distortion_irrep_numbers=distortion_irrep_numbers,
        include_flags=include_flags,
        mode_counts=mode_counts,
        coefficient_summaries=coefficient_summaries,
        mode_summaries=mode_summaries,
        primary_irreps=primary_irreps,
        modes_irrep_numbers=modes_irrep_numbers,
        irrep_count=int(sum(_numbers(modes.get("irrepCount")) or _numbers(distortion.get("irrepCount")) or [0])),
        modulation_count=modulation_count,
        subgroup_cell_size=subgroup_cell_size,
        classification=_classification(include_flags, modulation_count, len(distortion_irrep_numbers)),
        associated_files=associated_files,
        tags={
            "modesFileName": (distortion.get("modesFileName") or "").strip(),
            "isoFileName": (distortion.get("isoFileName") or "").strip(),
            "atomsFileName": (distortion.get("atomsFileName") or "").strip(),
        },
        child_basic_structure=child_basic_structure,
        internal_checks={
            "associated_files_all_present": not missing_associated_files,
            "missing_associated_files": missing_associated_files,
            "coefficient_lengths_consistent": coefficient_consistency_ok,
            "primary_irrep_alignment": primary_irrep_alignment,
            "atoms_section_present": child_basic_structure is not None,
        },
    )

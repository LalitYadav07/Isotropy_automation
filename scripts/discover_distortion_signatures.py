from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import time
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

try:
    import matplotlib.pyplot as plt
    import numpy as np
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.core import Lattice, Structure
    from pymatgen.io.cif import CifParser, CifWriter
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
except ImportError as exc:
    raise SystemExit(
        "Missing dependency. Run this script with the project virtualenv, for example:\n"
        "  ./.venv/bin/python isotropy_project/scripts/discover_distortion_signatures.py <parent.cif>"
    ) from exc

try:
    from _isotropy_paths import DISCOVERY_RUNS_DIR, configure_environment, resolve_input_path
except ModuleNotFoundError:
    from isotropy_project.scripts._isotropy_paths import DISCOVERY_RUNS_DIR, configure_environment, resolve_input_path


configure_environment()


PAGE_LENGTH = 1000
SCREEN_WIDTH = 240
ISO_TIMEOUT_SECONDS = 5
COUPLED_TIMEOUT_SECONDS = 30
PARENT_SYMPREC = 1e-3
VERIFY_SYMPREC = 1e-2
PROBE_AMPLITUDE = 0.02
XRD_WAVELENGTH = "CuKa"
REPORT_2THETA_MIN = 10
REPORT_2THETA_MAX = 90
PEAK_MATCH_TOLERANCE = 0.18
EMERGENT_PEAK_INTENSITY_THRESHOLD = 0.2
PARENT_PEAK_INTENSITY_THRESHOLD = 1.0
CHANGED_PEAK_DELTA_THRESHOLD = 0.5
MODELING_DIRECTION_PREFIX = "P"
EXACT_MATCH_DISTANCE_TOL = 5e-4
GOOD_MATCH_DISTANCE_TOL = 5e-3
POOR_MATCH_DISTANCE_TOL = 5e-2


@dataclass
class SiteSummary:
    index: int
    species: str
    wyckoff: str
    frac_coords: list[float]


@dataclass
class ParentInfo:
    source_cif: str
    standardized_cif: str
    formula: str
    sg_num: int
    sg_symbol: str
    lattice_parameters: list[float]
    site_summaries: list[SiteSummary]
    occupied_wyckoffs: list[str]
    raw_sg_num: int
    raw_sg_symbol: str


@dataclass
class KPointInfo:
    label: str
    coordinates: str
    is_fixed: bool


@dataclass
class BranchInfo:
    branch_id: str
    subgroup_key: str
    kpoint: str
    kpoint_coordinates: str
    irrep: str
    direction: str
    order_parameter: str
    sg_num: int
    symbol: str
    basis_text: str
    basis_matrix: list[list[float]]
    origin_text: str
    origin_vector: list[float]
    is_primary_direction: bool


class IsoSession:
    def __init__(self) -> None:
        self.call_count = 0

    def run(self, commands: str, *, timeout_seconds: int = ISO_TIMEOUT_SECONDS, prompt_returns: int = 0) -> str:
        self.call_count += 1
        full_cmd = (
            f"PAGE {PAGE_LENGTH}\n"
            f"SCREEN {SCREEN_WIDTH}\n"
            f"{commands.strip()}\n"
            f"{chr(10) * prompt_returns}"
            "QUIT\n"
        )
        try:
            result = subprocess.run(
                ["iso"],
                input=full_cmd,
                capture_output=True,
                text=True,
                check=False,
                timeout=timeout_seconds,
            )
        except subprocess.TimeoutExpired as exc:
            return exc.stdout or ""
        return result.stdout


def parse_number(token: str) -> float:
    cleaned = token.strip().replace("−", "-")
    try:
        return float(Fraction(cleaned))
    except (ValueError, ZeroDivisionError):
        return float(cleaned)


def parse_vector(text: str) -> list[float]:
    return [parse_number(part) for part in text.strip().strip("()").split(",")]


def parse_basis_vectors(text: str) -> list[list[float]]:
    return [parse_vector(match) for match in re.findall(r"\(([^)]+)\)", text)]


def cell_command(basis_matrix: list[list[float]]) -> str:
    rows = []
    for row in basis_matrix:
        rows.append(",".join(str(Fraction(value).limit_denominator()) for value in row))
    return " ".join(rows)


def periodic_distance(frac_a: np.ndarray, frac_b: np.ndarray) -> float:
    delta = np.mod(frac_a - frac_b + 0.5, 1.0) - 0.5
    return float(np.linalg.norm(delta))


def transform_parent_point_to_child(
    point_parent: list[float],
    basis_matrix: list[list[float]],
    origin_vector: list[float],
) -> np.ndarray:
    basis = np.array(basis_matrix, dtype=float)
    inverse = np.linalg.inv(basis)
    origin = np.array(origin_vector, dtype=float)
    return np.matmul(np.array(point_parent, dtype=float) - origin, inverse)


def transform_parent_vector_to_child(vector_parent: list[float], basis_matrix: list[list[float]]) -> np.ndarray:
    basis = np.array(basis_matrix, dtype=float)
    inverse = np.linalg.inv(basis)
    return np.matmul(np.array(vector_parent, dtype=float), inverse)


def structure_signature(structure: Structure) -> str:
    analyzer = SpacegroupAnalyzer(structure, symprec=VERIFY_SYMPREC)
    symmetrized = analyzer.get_symmetrized_structure()
    positions = []
    for site in symmetrized:
        coords = ",".join(f"{value:.5f}" for value in np.mod(site.frac_coords, 1.0))
        positions.append(f"{site.species_string}:{coords}")
    positions.sort()
    return f"{analyzer.get_space_group_number()}|{analyzer.get_space_group_symbol()}|{'|'.join(positions)}"


def standardize_parent(parent_cif: Path, run_dir: Path) -> tuple[ParentInfo, Structure]:
    raw_structure = CifParser(str(parent_cif)).parse_structures(primitive=False)[0]
    raw_analyzer = SpacegroupAnalyzer(raw_structure, symprec=PARENT_SYMPREC)

    standard_structure = raw_analyzer.get_conventional_standard_structure()
    standard_path = run_dir / "parent_standardized.cif"
    CifWriter(standard_structure).write_file(standard_path)

    analyzer = SpacegroupAnalyzer(standard_structure, symprec=PARENT_SYMPREC)
    dataset = analyzer.get_symmetry_dataset()

    site_summaries: list[SiteSummary] = []
    occupied_wyckoffs: list[str] = []
    for index, (site, wyckoff_letter) in enumerate(zip(standard_structure, dataset.wyckoffs, strict=False)):
        wyckoff = str(wyckoff_letter).lower()
        occupied_wyckoffs.append(wyckoff)
        site_summaries.append(
            SiteSummary(
                index=index,
                species=site.species_string,
                wyckoff=wyckoff,
                frac_coords=[float(value) for value in np.mod(site.frac_coords, 1.0)],
            )
        )

    parent_info = ParentInfo(
        source_cif=str(parent_cif),
        standardized_cif=str(standard_path),
        formula=standard_structure.composition.reduced_formula,
        sg_num=analyzer.get_space_group_number(),
        sg_symbol=analyzer.get_space_group_symbol(),
        lattice_parameters=[float(value) for value in standard_structure.lattice.abc + standard_structure.lattice.angles],
        site_summaries=site_summaries,
        occupied_wyckoffs=sorted(set(occupied_wyckoffs)),
        raw_sg_num=raw_analyzer.get_space_group_number(),
        raw_sg_symbol=raw_analyzer.get_space_group_symbol(),
    )
    return parent_info, standard_structure


def get_kpoints(iso: IsoSession, parent_sg: int) -> list[KPointInfo]:
    out = iso.run(
        f"""
        VALUE PARENT {parent_sg}
        SHOW KPOINT
        DISPLAY KPOINT
        """
    )
    kpoints: list[KPointInfo] = []
    for line in out.splitlines():
        match = re.match(r"^([A-Z0-9]+)\s+\(([^)]+)\)", line.strip())
        if not match:
            continue
        coordinates = match.group(2)
        kpoints.append(
            KPointInfo(
                label=match.group(1),
                coordinates=coordinates,
                is_fixed=re.search(r"[a-z]", coordinates) is None,
            )
        )
    return kpoints


def get_irreps(iso: IsoSession, parent_sg: int, kpoint: str) -> list[str]:
    out = iso.run(
        f"""
        VALUE PARENT {parent_sg}
        VALUE KPOINT {kpoint}
        SHOW IRREP
        DISPLAY IRREP
        """
    )
    irreps: list[str] = []
    header_seen = False
    for line in out.splitlines():
        stripped = line.strip()
        if "Irrep" in stripped:
            header_seen = True
            continue
        if not header_seen or not stripped or stripped == "*":
            continue
        token = stripped.split()[0]
        if token != "*":
            irreps.append(token)
    return irreps


def get_subgroups(iso: IsoSession, parent_sg: int, kpoint: str, irrep: str) -> list[dict[str, object]]:
    out = iso.run(
        f"""
        VALUE PARENT {parent_sg}
        VALUE KPOINT {kpoint}
        VALUE IRREP {irrep}
        SHOW SUBGROUP
        SHOW DIRECTION
        SHOW DIRECTION VECTOR
        SHOW BASIS
        SHOW ORIGIN
        DISPLAY ISOTROPY
        """
    )
    pattern = re.compile(
        r"^(\d+)\s+(\S+)\s+(\S+)\s+(\(.+?\))\s+(\([^)]*\),\([^)]*\),\([^)]*\))\s+(\([^)]*\))$"
    )
    subgroups: list[dict[str, object]] = []
    for line in out.splitlines():
        match = pattern.match(line.strip())
        if not match:
            continue
        subgroups.append(
            {
                "sg_num": int(match.group(1)),
                "symbol": match.group(2),
                "direction": match.group(3),
                "order_parameter": match.group(4),
                "basis_text": match.group(5),
                "basis_matrix": parse_basis_vectors(match.group(5)),
                "origin_text": match.group(6),
                "origin_vector": parse_vector(match.group(6)),
            }
        )
    return subgroups


def get_domain_count(iso: IsoSession, parent_sg: int, kpoint: str, irrep: str, direction: str) -> int | None:
    out = iso.run(
        f"""
        VALUE PARENT {parent_sg}
        VALUE KPOINT {kpoint}
        VALUE IRREP {irrep}
        VALUE DIRECTION {direction}
        SHOW DOMAINS
        DISPLAY ISOTROPY
        """
    )
    domain_lines = [line.strip() for line in out.splitlines() if line.strip().isdigit()]
    return len(domain_lines) if domain_lines else None


def get_vector_distortion_rows(
    iso: IsoSession,
    parent_sg: int,
    branch: BranchInfo,
    occupied_wyckoffs: list[str],
) -> list[tuple[str, list[float], list[float]]]:
    wyckoff_str = " ".join(letter.upper() for letter in occupied_wyckoffs)
    out = iso.run(
        f"""
        VALUE PARENT {parent_sg}
        VALUE KPOINT {branch.kpoint}
        VALUE IRREP {branch.irrep}
        VALUE DIRECTION {branch.direction}
        VALUE CELL {cell_command(branch.basis_matrix)}
        VALUE WYCKOFF {wyckoff_str}
        SHOW WYCKOFF
        SHOW MICROSCOPIC VECTOR
        DISPLAY DISTORTION
        """
    )
    rows: list[tuple[str, list[float], list[float]]] = []
    current_wyckoff: str | None = None
    pattern = re.compile(r"^([A-Za-z0-9]+)?\s*\(([^)]+)\)\s+\(([^)]+)\)")
    for line in out.splitlines():
        match = pattern.match(line.rstrip())
        if not match:
            continue
        if match.group(1):
            current_wyckoff = match.group(1).lower()
        if current_wyckoff is None:
            continue
        rows.append(
            (
                current_wyckoff,
                [parse_number(part) for part in match.group(2).split(",")],
                [parse_number(part) for part in match.group(3).split(",")],
            )
        )
    return rows


def screen_irrep_channels(
    iso: IsoSession,
    parent_info: ParentInfo,
    branches: list[BranchInfo],
) -> tuple[list[dict[str, object]], dict[str, list[tuple[str, list[float], list[float]]]]]:
    grouped: dict[tuple[str, str], list[BranchInfo]] = {}
    for branch in branches:
        if branch.is_primary_direction:
            grouped.setdefault((branch.kpoint, branch.irrep), []).append(branch)

    rows_cache: dict[str, list[tuple[str, list[float], list[float]]]] = {}
    summaries: list[dict[str, object]] = []
    for (kpoint, irrep), group in sorted(grouped.items()):
        ordered = sorted(group, key=branch_sort_key)
        vector_branch = None
        sample_row_count = 0
        attempted_branch_ids = []
        for branch in ordered:
            rows = get_vector_distortion_rows(iso, parent_info.sg_num, branch, parent_info.occupied_wyckoffs)
            rows_cache[branch.branch_id] = rows
            attempted_branch_ids.append(branch.branch_id)
            if rows:
                vector_branch = branch
                sample_row_count = len(rows)
                break
        summaries.append(
            {
                "screen_id": f"I{len(summaries) + 1:03d}",
                "kpoint": kpoint,
                "irrep": irrep,
                "kpoint_coordinates": ordered[0].kpoint_coordinates,
                "primary_branch_count": len(ordered),
                "attempted_branch_ids": attempted_branch_ids,
                "has_vector_channel": vector_branch is not None,
                "screening_status": "vector_active" if vector_branch is not None else "no_vector_channel",
                "exemplar_branch_id": vector_branch.branch_id if vector_branch is not None else None,
                "exemplar_direction": vector_branch.direction if vector_branch is not None else None,
                "sample_row_count": sample_row_count,
            }
        )
    return summaries, rows_cache


def get_coupled_subgroups(
    iso: IsoSession,
    parent_sg: int,
    irrep_a: str,
    irrep_b: str,
) -> list[dict[str, object]]:
    out = iso.run(
        f"""
        VALUE PARENT {parent_sg}
        VALUE IRREP {irrep_a} {irrep_b}
        SHOW SUBGROUP
        SHOW DIRECTION
        SHOW DIRECTION VECTOR
        SHOW BASIS
        SHOW ORIGIN
        SHOW IRREP
        DISPLAY ISOTROPY COUPLED
        """,
        timeout_seconds=COUPLED_TIMEOUT_SECONDS,
        prompt_returns=1,
    )
    pattern = re.compile(
        r"^(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\(.+?\))\s+(\([^)]*\),\([^)]*\),\([^)]*\))\s+(\([^)]*\))$"
    )
    subgroups: list[dict[str, object]] = []
    for line in out.splitlines():
        match = pattern.match(line.strip())
        if not match:
            continue
        subgroups.append(
            {
                "reducible_irrep": match.group(1),
                "sg_num": int(match.group(2)),
                "symbol": match.group(3),
                "direction": match.group(4),
                "order_parameter": match.group(5),
                "basis_text": match.group(6),
                "basis_matrix": parse_basis_vectors(match.group(6)),
                "origin_text": match.group(7),
                "origin_vector": parse_vector(match.group(7)),
            }
        )
    return subgroups


def transition_class_rows(report_settings: dict[str, object], stats: dict[str, object]) -> list[dict[str, str]]:
    coupled_status = "cataloged" if report_settings["include_coupled_catalog"] else "not_run"
    coupled_detail = (
        (
            f"{stats['coupled_subgroup_count']} coupled subgroup branches from {stats['coupled_pair_count']} expanded irrep pairs "
            f"out of {stats['possible_coupled_pair_count']} possible screened vector-active pairs."
        )
        if report_settings["include_coupled_catalog"]
        else "Coupled-primary intersections were not requested in this run."
    )
    return [
        {
            "class": "Single-irrep commensurate fixed-k displacive",
            "status": "implemented",
            "detail": "Enumerated rigorously and probed with canonical vector distortions.",
        },
        {
            "class": "Coupled commensurate fixed-k displacive",
            "status": coupled_status,
            "detail": coupled_detail,
        },
        {
            "class": "Parameterized commensurate k-lines / k-planes",
            "status": "documented_only",
            "detail": f"{stats['parameterized_kpoint_count']} parameterized k manifolds are listed but not expanded automatically.",
        },
        {
            "class": "One-arm and kernel-only searches",
            "status": "not_yet_automated",
            "detail": "Supported by ISOTROPY commands, but not yet scripted in this workflow.",
        },
        {
            "class": "Incommensurate superspace transitions",
            "status": "not_yet_automated",
            "detail": "These require superspace treatment rather than the current commensurate child-cell probe builder.",
        },
        {
            "class": "Occupational / scalar structural order",
            "status": "not_yet_automated",
            "detail": "Current modeling uses microscopic vectors only; scalar order parameters are catalog-only at present.",
        },
        {
            "class": "Reconstructive transitions",
            "status": "requires_child_structure",
            "detail": "Need both phases and a common-subgroup search with comsubs; not discoverable from a parent CIF alone.",
        },
    ]


def subgroup_key(data: dict[str, object]) -> str:
    return f"{data['sg_num']}|{data['symbol']}|{data['basis_text']}|{data['origin_text']}"


def branch_sort_key(branch: BranchInfo) -> tuple[int, str, str, str]:
    return (0 if branch.is_primary_direction else 1, branch.kpoint, branch.irrep, branch.direction)


def build_branch_catalog(
    iso: IsoSession,
    parent_info: ParentInfo,
    kpoints: list[KPointInfo],
) -> tuple[list[BranchInfo], list[KPointInfo]]:
    fixed_kpoints = [kpoint for kpoint in kpoints if kpoint.is_fixed]
    branches: list[BranchInfo] = []
    counter = 1
    for kpoint in fixed_kpoints:
        for irrep in get_irreps(iso, parent_info.sg_num, kpoint.label):
            for subgroup in get_subgroups(iso, parent_info.sg_num, kpoint.label, irrep):
                key = subgroup_key(subgroup)
                branches.append(
                    BranchInfo(
                        branch_id=f"B{counter:04d}",
                        subgroup_key=key,
                        kpoint=kpoint.label,
                        kpoint_coordinates=kpoint.coordinates,
                        irrep=irrep,
                        direction=str(subgroup["direction"]),
                        order_parameter=str(subgroup["order_parameter"]),
                        sg_num=int(subgroup["sg_num"]),
                        symbol=str(subgroup["symbol"]),
                        basis_text=str(subgroup["basis_text"]),
                        basis_matrix=[[float(value) for value in row] for row in subgroup["basis_matrix"]],
                        origin_text=str(subgroup["origin_text"]),
                        origin_vector=[float(value) for value in subgroup["origin_vector"]],
                        is_primary_direction=str(subgroup["direction"]).startswith(MODELING_DIRECTION_PREFIX),
                    )
                )
                counter += 1
    return branches, fixed_kpoints


def basis_det(basis_matrix: list[list[float]]) -> float:
    return float(np.linalg.det(np.array(basis_matrix, dtype=float)))


def build_child_structure(
    parent_structure: Structure,
    parent_wyckoffs: list[str],
    basis_matrix: list[list[float]],
    origin_vector: list[float],
    max_extent: int = 5,
    tolerance: float = 1e-5,
) -> tuple[Structure, dict[str, object], list[dict[str, object]]]:
    basis = np.array(basis_matrix, dtype=float)
    origin = np.array(origin_vector, dtype=float)
    inverse = np.linalg.inv(basis)
    child_lattice = Lattice(np.matmul(basis, parent_structure.lattice.matrix))
    expected_site_count = max(1, int(round(abs(np.linalg.det(basis)) * len(parent_structure))))

    final_species: list[str] = []
    final_coords: list[list[float]] = []
    final_metadata: list[dict[str, object]] = []
    final_extent = None

    for extent in range(2, max_extent + 1):
        species: list[str] = []
        coords: list[list[float]] = []
        metadata: list[dict[str, object]] = []
        for parent_index, site in enumerate(parent_structure):
            parent_wyckoff = parent_wyckoffs[parent_index] if parent_index < len(parent_wyckoffs) else "?"
            for translation in product(range(-extent, extent + 1), repeat=3):
                parent_coords = np.array(site.frac_coords) + np.array(translation, dtype=float)
                child_coords = np.matmul(parent_coords - origin, inverse)
                if np.all(child_coords >= -tolerance) and np.all(child_coords < 1 - tolerance):
                    wrapped = np.mod(child_coords, 1.0)
                    is_duplicate = False
                    for existing_species, existing_coords in zip(species, coords, strict=False):
                        if existing_species != site.species_string:
                            continue
                        if periodic_distance(np.array(existing_coords), wrapped) < 5e-4:
                            is_duplicate = True
                            break
                    if not is_duplicate:
                        species.append(site.species_string)
                        coords.append([float(value) for value in wrapped])
                        metadata.append(
                            {
                                "species": site.species_string,
                                "parent_index": parent_index,
                                "wyckoff": parent_wyckoff,
                            }
                        )
        if len(species) == expected_site_count:
            final_species = species
            final_coords = coords
            final_metadata = metadata
            final_extent = extent
            break
        if len(species) > len(final_species):
            final_species = species
            final_coords = coords
            final_metadata = metadata
            final_extent = extent

    structure = Structure(child_lattice, final_species, final_coords, coords_are_cartesian=False)
    return structure, {
        "expected_site_count": expected_site_count,
        "actual_site_count": len(final_species),
        "translation_extent": final_extent,
        "determinant": basis_det(basis_matrix),
    }, final_metadata


def aggregate_distortion_rows(
    rows: list[tuple[str, list[float], list[float]]],
    basis_matrix: list[list[float]],
    origin_vector: list[float],
) -> list[dict[str, object]]:
    aggregated: list[dict[str, object]] = []
    for row_index, (wyckoff, point_parent, vector_parent) in enumerate(rows):
        point_child = np.mod(
            transform_parent_point_to_child(point_parent, basis_matrix, origin_vector),
            1.0,
        )
        vector_child = transform_parent_vector_to_child(vector_parent, basis_matrix)

        matched = None
        for entry in aggregated:
            if entry["wyckoff"] != wyckoff:
                continue
            if periodic_distance(np.array(entry["point_child"], dtype=float), point_child) <= EXACT_MATCH_DISTANCE_TOL:
                matched = entry
                break

        if matched is None:
            aggregated.append(
                {
                    "wyckoff": wyckoff,
                    "point_child": [float(value) for value in point_child],
                    "vector_child": [float(value) for value in vector_child],
                    "source_row_indices": [row_index],
                    "source_points_parent": [[float(value) for value in point_parent]],
                    "source_vectors_parent": [[float(value) for value in vector_parent]],
                }
            )
            continue

        matched["vector_child"] = [
            float(value)
            for value in (
                np.array(matched["vector_child"], dtype=float) + vector_child
            )
        ]
        matched["source_row_indices"].append(row_index)
        matched["source_points_parent"].append([float(value) for value in point_parent])
        matched["source_vectors_parent"].append([float(value) for value in vector_parent])

    return aggregated


def match_rows_to_structure(
    rows: list[dict[str, object]],
    structure: Structure,
    site_metadata: list[dict[str, object]],
) -> dict[str, object]:
    used_indices: set[int] = set()
    assignments = []
    unmatched_rows = 0
    max_distance = 0.0
    wyckoff_fallback_count = 0
    for row_index, row in enumerate(rows):
        point_array = np.mod(np.array(row["point_child"], dtype=float), 1.0)
        best_match: tuple[float, int] | None = None
        preferred_indices = [
            site_index
            for site_index, metadata in enumerate(site_metadata)
            if site_index not in used_indices and metadata["wyckoff"] == row["wyckoff"]
        ]
        fallback_used = False
        if not preferred_indices:
            preferred_indices = [site_index for site_index in range(len(structure)) if site_index not in used_indices]
            fallback_used = True
        for site_index in preferred_indices:
            site = structure[site_index]
            distance = periodic_distance(np.mod(site.frac_coords, 1.0), point_array)
            if best_match is None or distance < best_match[0]:
                best_match = (distance, site_index)
        if best_match is None:
            unmatched_rows += 1
            continue
        used_indices.add(best_match[1])
        max_distance = max(max_distance, best_match[0])
        if fallback_used:
            wyckoff_fallback_count += 1
        assignments.append(
            {
                "row_index": row_index,
                "site_index": best_match[1],
                "distance": best_match[0],
                "point": [float(value) for value in point_array],
                "vector": [float(value) for value in row["vector_child"]],
                "wyckoff": row["wyckoff"],
                "source_row_indices": row["source_row_indices"],
            }
        )
    return {
        "matched_count": len(assignments),
        "unmatched_count": unmatched_rows,
        "max_distance": max_distance,
        "wyckoff_fallback_count": wyckoff_fallback_count,
        "assignments": assignments,
    }


def select_reference_structure(
    parent_structure: Structure,
    parent_wyckoffs: list[str],
    branch: BranchInfo,
    rows: list[tuple[str, list[float], list[float]]],
) -> tuple[Structure, dict[str, object], list[dict[str, object]], dict[str, object], list[dict[str, object]], str]:
    candidates = []
    for origin_mode, origin_vector in (
        ("subgroup_origin", branch.origin_vector),
        ("zero_origin", [0.0, 0.0, 0.0]),
    ):
        structure, build_info, site_metadata = build_child_structure(
            parent_structure,
            parent_wyckoffs,
            branch.basis_matrix,
            origin_vector,
        )
        transformed_rows = aggregate_distortion_rows(rows, branch.basis_matrix, origin_vector)
        match_info = match_rows_to_structure(transformed_rows, structure, site_metadata)
        candidates.append((structure, build_info, site_metadata, match_info, transformed_rows, origin_mode))

    candidates.sort(
        key=lambda item: (
            item[3]["unmatched_count"],
            item[3]["wyckoff_fallback_count"],
            item[3]["max_distance"],
            abs(item[1]["expected_site_count"] - item[1]["actual_site_count"]),
            -item[3]["matched_count"],
        )
    )
    return candidates[0]


def apply_probe_distortion(
    structure: Structure,
    match_info: dict[str, object],
    rows: list[dict[str, object]],
    amplitude: float,
) -> tuple[Structure, dict[str, object]]:
    distorted = structure.copy()
    species = []
    norms = []
    moved_sites = 0
    for assignment in match_info["assignments"]:
        site_index = assignment["site_index"]
        vector = np.array(rows[assignment["row_index"]]["vector_child"], dtype=float)
        if np.linalg.norm(vector) <= 1e-12:
            continue
        moved_sites += 1
        norms.append(float(np.linalg.norm(vector)))
        species.append(distorted[site_index].species_string)
        new_coords = np.mod(distorted[site_index].frac_coords + amplitude * vector, 1.0)
        distorted.replace(site_index, distorted[site_index].species, new_coords, coords_are_cartesian=False)

    return distorted, {
        "moved_site_count": moved_sites,
        "moved_species": sorted(set(species)),
        "max_mode_norm": max(norms, default=0.0),
        "rms_mode_norm": float(np.sqrt(np.mean(np.square(norms)))) if norms else 0.0,
    }


def verify_structure(structure: Structure) -> tuple[int, str]:
    analyzer = SpacegroupAnalyzer(structure, symprec=VERIFY_SYMPREC)
    return analyzer.get_space_group_number(), analyzer.get_space_group_symbol()


def strongest_peaks(pattern, count: int = 12) -> list[dict[str, float]]:
    peaks = [
        (two_theta, intensity)
        for two_theta, intensity in zip(pattern.x, pattern.y, strict=False)
        if REPORT_2THETA_MIN <= two_theta <= REPORT_2THETA_MAX
    ]
    peaks.sort(key=lambda item: item[1], reverse=True)
    selected = []
    for two_theta, intensity in peaks:
        selected.append({"two_theta": float(two_theta), "intensity": float(intensity)})
        if len(selected) == count:
            break
    return selected


def filtered_peaks(pattern, min_intensity: float) -> list[dict[str, float]]:
    peaks = []
    for two_theta, intensity in zip(pattern.x, pattern.y, strict=False):
        if REPORT_2THETA_MIN <= two_theta <= REPORT_2THETA_MAX and intensity >= min_intensity:
            peaks.append({"two_theta": float(two_theta), "intensity": float(intensity)})
    return peaks


def emergent_peaks(parent_pattern, distorted_pattern) -> list[dict[str, float]]:
    parent_peaks = filtered_peaks(parent_pattern, min_intensity=EMERGENT_PEAK_INTENSITY_THRESHOLD)
    distorted_peaks = filtered_peaks(distorted_pattern, min_intensity=EMERGENT_PEAK_INTENSITY_THRESHOLD)
    results = []
    for peak in distorted_peaks:
        if not any(abs(peak["two_theta"] - ref["two_theta"]) <= PEAK_MATCH_TOLERANCE for ref in parent_peaks):
            results.append(peak)
    results.sort(key=lambda item: (-item["intensity"], item["two_theta"]))
    return sorted(results[:12], key=lambda item: item["two_theta"])


def suppressed_parent_peaks(parent_pattern, distorted_pattern) -> list[dict[str, float]]:
    parent_peaks = filtered_peaks(parent_pattern, min_intensity=PARENT_PEAK_INTENSITY_THRESHOLD)
    distorted_peaks = filtered_peaks(distorted_pattern, min_intensity=EMERGENT_PEAK_INTENSITY_THRESHOLD)
    results = []
    for peak in parent_peaks:
        if not any(abs(peak["two_theta"] - ref["two_theta"]) <= PEAK_MATCH_TOLERANCE for ref in distorted_peaks):
            results.append(peak)
    results.sort(key=lambda item: (-item["intensity"], item["two_theta"]))
    return sorted(results[:12], key=lambda item: item["two_theta"])


def changed_parent_peaks(parent_pattern, distorted_pattern) -> list[dict[str, float]]:
    parent_peaks = filtered_peaks(parent_pattern, min_intensity=PARENT_PEAK_INTENSITY_THRESHOLD)
    distorted_peaks = filtered_peaks(distorted_pattern, min_intensity=EMERGENT_PEAK_INTENSITY_THRESHOLD)
    results = []
    for peak in parent_peaks:
        matches = [ref for ref in distorted_peaks if abs(peak["two_theta"] - ref["two_theta"]) <= PEAK_MATCH_TOLERANCE]
        if not matches:
            continue
        closest = min(matches, key=lambda ref: abs(peak["two_theta"] - ref["two_theta"]))
        delta = closest["intensity"] - peak["intensity"]
        if abs(delta) >= CHANGED_PEAK_DELTA_THRESHOLD:
            results.append(
                {
                    "two_theta": peak["two_theta"],
                    "parent_intensity": peak["intensity"],
                    "distorted_intensity": closest["intensity"],
                    "delta_intensity": delta,
                }
            )
    results.sort(key=lambda item: abs(item["delta_intensity"]), reverse=True)
    return results[:8]


def signature_key(signature: dict[str, object]) -> str:
    emergent = tuple(round(item["two_theta"], 1) for item in signature["emergent_peaks"][:4])
    changed = tuple(
        f"{round(item['two_theta'], 1)}:{'+' if item['delta_intensity'] >= 0 else '-'}"
        for item in signature["changed_parent_peaks"][:3]
    )
    return json.dumps(
        {
            "emergent": emergent,
            "changed": changed,
            "suppressed": tuple(round(item["two_theta"], 1) for item in signature["suppressed_parent_peaks"][:3]),
        },
        sort_keys=True,
    )


def make_signature_plot(
    parent_structure: Structure,
    distorted_structure: Structure,
    label: str,
    plot_dir: Path,
) -> tuple[str, dict[str, object]]:
    calculator = XRDCalculator(wavelength=XRD_WAVELENGTH)
    parent_pattern = calculator.get_pattern(parent_structure)
    distorted_pattern = calculator.get_pattern(distorted_structure)

    plt.figure(figsize=(10, 6))
    plt.vlines(parent_pattern.x, 0, parent_pattern.y, color="#456990", alpha=0.7, label="Parent")
    plt.vlines(distorted_pattern.x, 0, distorted_pattern.y, color="#ef8354", alpha=0.7, label="Probe distortion")
    plt.xlabel("2θ (degrees)")
    plt.ylabel("Relative intensity")
    plt.xlim(REPORT_2THETA_MIN, REPORT_2THETA_MAX)
    plt.title(label)
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()

    filename = re.sub(r"[^A-Za-z0-9_.-]+", "_", label) + ".png"
    path = plot_dir / filename
    plt.savefig(path, dpi=180)
    plt.close()

    summary = {
        "parent_strongest": strongest_peaks(parent_pattern),
        "probe_strongest": strongest_peaks(distorted_pattern),
        "emergent_peaks": emergent_peaks(parent_pattern, distorted_pattern),
        "suppressed_parent_peaks": suppressed_parent_peaks(parent_pattern, distorted_pattern),
        "changed_parent_peaks": changed_parent_peaks(parent_pattern, distorted_pattern),
    }
    summary["signature_key"] = signature_key(summary)
    return filename, summary


def write_structure(structure: Structure, name: str, structure_dir: Path) -> str:
    filename = re.sub(r"[^A-Za-z0-9_.-]+", "_", name) + ".cif"
    path = structure_dir / filename
    CifWriter(structure).write_file(path)
    return filename


def parse_residual_peaks(path: Path) -> list[dict[str, float]]:
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        peaks = []
        for entry in data:
            if isinstance(entry, dict):
                two_theta = entry.get("two_theta") or entry.get("angle") or entry.get("x")
                intensity = entry.get("intensity") or entry.get("y") or 1.0
            else:
                two_theta, intensity = entry[0], entry[1] if len(entry) > 1 else 1.0
            peaks.append({"two_theta": float(two_theta), "intensity": float(intensity)})
        return peaks

    peaks = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue
            try:
                peaks.append({"two_theta": float(row[0]), "intensity": float(row[1]) if len(row) > 1 else 1.0})
            except ValueError:
                continue
    return peaks


def residual_match_score(signature: dict[str, object], residual_peaks: list[dict[str, float]]) -> dict[str, object]:
    predicted = signature["emergent_peaks"] + [
        {"two_theta": item["two_theta"], "intensity": abs(item["delta_intensity"])}
        for item in signature["changed_parent_peaks"]
    ]
    if not residual_peaks or not predicted:
        return {"score": 0.0, "matched_residual_count": 0, "matched_predicted_count": 0}

    matched_residual = set()
    matched_predicted = set()
    for predicted_index, predicted_peak in enumerate(predicted):
        for residual_index, residual_peak in enumerate(residual_peaks):
            if abs(predicted_peak["two_theta"] - residual_peak["two_theta"]) <= PEAK_MATCH_TOLERANCE:
                matched_residual.add(residual_index)
                matched_predicted.add(predicted_index)
    coverage = len(matched_residual) / len(residual_peaks)
    precision = len(matched_predicted) / len(predicted)
    return {
        "score": round(0.7 * coverage + 0.3 * precision, 4),
        "matched_residual_count": len(matched_residual),
        "matched_predicted_count": len(matched_predicted),
        "coverage": round(coverage, 4),
        "precision": round(precision, 4),
    }


def embedding_quality(match_info: dict[str, object], build_info: dict[str, object]) -> str:
    site_count_error = abs(build_info["actual_site_count"] - build_info["expected_site_count"])
    if match_info["unmatched_count"] > 0 or site_count_error > 0:
        return "failed"
    if match_info["wyckoff_fallback_count"] > 0:
        return "partial"
    if match_info["max_distance"] <= EXACT_MATCH_DISTANCE_TOL:
        return "exact"
    if match_info["max_distance"] <= GOOD_MATCH_DISTANCE_TOL:
        return "good"
    if match_info["max_distance"] <= POOR_MATCH_DISTANCE_TOL:
        return "poor"
    return "failed"


def unique_subgroup_summary(branches: list[BranchInfo], domain_counts: dict[str, int | None]) -> list[dict[str, object]]:
    grouped: dict[str, list[BranchInfo]] = {}
    for branch in branches:
        grouped.setdefault(branch.subgroup_key, []).append(branch)

    rows = []
    for subgroup_id, (key, group) in enumerate(sorted(grouped.items()), start=1):
        ordered = sorted(group, key=branch_sort_key)
        exemplar = ordered[0]
        rows.append(
            {
                "subgroup_id": f"U{subgroup_id:03d}",
                "subgroup_key": key,
                "sg_num": exemplar.sg_num,
                "symbol": exemplar.symbol,
                "basis_text": exemplar.basis_text,
                "origin_text": exemplar.origin_text,
                "branch_count": len(group),
                "primary_branch_count": sum(1 for branch in group if branch.is_primary_direction),
                "example_branches": [branch.branch_id for branch in ordered[:6]],
                "domain_count": domain_counts.get(ordered[0].branch_id),
            }
        )
    return rows


def run_pipeline(args: argparse.Namespace) -> dict[str, object]:
    start_total = time.perf_counter()
    run_dir = DISCOVERY_RUNS_DIR / args.label
    structure_dir = run_dir / "candidate_structures"
    plot_dir = run_dir / "signature_plots"
    run_dir.mkdir(parents=True, exist_ok=True)
    structure_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    timings: dict[str, float] = {}

    t0 = time.perf_counter()
    parent_path = resolve_input_path(args.parent_cif)
    parent_info, parent_structure = standardize_parent(parent_path, run_dir)
    timings["parent_standardization_seconds"] = round(time.perf_counter() - t0, 3)

    iso = IsoSession()

    t0 = time.perf_counter()
    kpoints = get_kpoints(iso, parent_info.sg_num)
    branches, fixed_kpoints = build_branch_catalog(iso, parent_info, kpoints)
    timings["catalog_enumeration_seconds"] = round(time.perf_counter() - t0, 3)

    t0 = time.perf_counter()
    irrep_screen, rows_cache = screen_irrep_channels(iso, parent_info, branches)
    timings["irrep_channel_screening_seconds"] = round(time.perf_counter() - t0, 3)

    domain_counts: dict[str, int | None] = {}
    grouped_branches = {}
    for branch in branches:
        grouped_branches.setdefault(branch.subgroup_key, []).append(branch)

    unique_subgroups = unique_subgroup_summary(branches, domain_counts)
    subgroup_lookup = {row["subgroup_key"]: row["subgroup_id"] for row in unique_subgroups}

    residual_peaks = parse_residual_peaks(resolve_input_path(args.residual_peaks)) if args.residual_peaks else []

    modeled_candidates: list[dict[str, object]] = []
    signature_groups: dict[str, list[str]] = {}
    seen_structures: dict[str, str] = {}

    t0 = time.perf_counter()
    primary_branches = [branch for branch in branches if branch.is_primary_direction]
    if args.max_modeled:
        primary_branches = primary_branches[: args.max_modeled]

    for candidate_index, branch in enumerate(primary_branches, start=1):
        candidate_id = f"C{candidate_index:03d}"
        notes: list[str] = []
        branch_domain_count = get_domain_count(iso, parent_info.sg_num, branch.kpoint, branch.irrep, branch.direction)
        domain_counts[branch.branch_id] = branch_domain_count

        rows = rows_cache.get(branch.branch_id)
        if rows is None:
            rows = get_vector_distortion_rows(iso, parent_info.sg_num, branch, parent_info.occupied_wyckoffs)
            rows_cache[branch.branch_id] = rows
        if not rows:
            modeled_candidates.append(
                {
                    "candidate_id": candidate_id,
                    "branch_id": branch.branch_id,
                    "subgroup_id": subgroup_lookup[branch.subgroup_key],
                    "status": "no_vector_mode",
                    "branch": asdict(branch),
                    "domain_count": branch_domain_count,
                    "embedding_quality": "not_applicable",
                    "notes": ["ISO did not return a microscopic vector distortion for the occupied parent Wyckoff sites."],
                }
            )
            continue

        reference_structure, build_info, site_metadata, match_info, transformed_rows, origin_mode = select_reference_structure(
            parent_structure,
            [site.wyckoff for site in parent_info.site_summaries],
            branch,
            rows,
        )
        unique_point_count = len(transformed_rows)
        if len(rows) != unique_point_count:
            notes.append(
                f"Collapsed {len(rows)} raw distortion rows into {unique_point_count} unique child-cell sites by summing repeated projected-vector contributions."
            )
        if match_info["matched_count"] != unique_point_count:
            notes.append(
                f"Only {match_info['matched_count']} of {unique_point_count} unique child-cell distortion sites were matched onto the constructed subgroup cell."
            )
        notes.append(f"Reference embedding used `{origin_mode}`.")
        if abs(build_info["actual_site_count"] - build_info["expected_site_count"]) > 0:
            notes.append(
                f"Constructed subgroup cell has {build_info['actual_site_count']} sites but expected approximately {build_info['expected_site_count']}."
            )
        if match_info["wyckoff_fallback_count"] > 0:
            notes.append(
                f"{match_info['wyckoff_fallback_count']} distortion-site assignments required a fallback beyond Wyckoff-label matching."
            )
        if match_info["max_distance"] > GOOD_MATCH_DISTANCE_TOL:
            notes.append(f"Maximum point-to-site matching distance was {match_info['max_distance']:.4f} in subgroup fractional coordinates.")

        distorted_structure, displacement_summary = apply_probe_distortion(
            reference_structure,
            match_info,
            transformed_rows,
            args.amplitude,
        )
        predicted_sg = f"{branch.sg_num} {branch.symbol}"
        verified_sg_num, verified_sg_symbol = verify_structure(distorted_structure)
        quality = embedding_quality(match_info, build_info)
        if verified_sg_num == branch.sg_num and quality in {"exact", "good"}:
            status = "verified"
        elif verified_sg_num == branch.sg_num and quality in {"partial", "poor"}:
            status = "tentative"
            notes.append(
                "The probe distortion reproduced the predicted subgroup, but the subgroup-cell embedding was not exact enough to treat this as a clean verification."
            )
        elif quality == "failed":
            status = "embedding_failed"
            notes.append(
                "The subgroup-cell embedding did not pass the geometric consistency checks, so this candidate should not be interpreted as a validated child structure."
            )
        else:
            status = "mismatch"
            notes.append(f"Canonical probe distortion verified as `{verified_sg_num} {verified_sg_symbol}` instead of predicted `{predicted_sg}`.")

        structure_filename = write_structure(distorted_structure, f"{candidate_id}_{branch.kpoint}_{branch.irrep}_{branch.direction}", structure_dir)
        plot_filename, signature = make_signature_plot(
            parent_structure,
            distorted_structure,
            f"{candidate_id} {branch.kpoint} {branch.irrep} {branch.direction}",
            plot_dir,
        )
        signature_groups.setdefault(signature["signature_key"], []).append(candidate_id)
        structure_sig = structure_signature(distorted_structure)
        duplicate_of = seen_structures.get(structure_sig)
        if duplicate_of is None:
            seen_structures[structure_sig] = candidate_id
        else:
            notes.append(f"Probe structure is symmetry-equivalent to candidate `{duplicate_of}` under the current signature metric.")

        residual_match = residual_match_score(signature, residual_peaks) if residual_peaks else None

        modeled_candidates.append(
            {
                "candidate_id": candidate_id,
                "branch_id": branch.branch_id,
                "subgroup_id": subgroup_lookup[branch.subgroup_key],
                "status": status,
                "domain_count": branch_domain_count,
                "branch": asdict(branch),
                "reference_build": build_info,
                "match": match_info,
                "raw_distortion_row_count": len(rows),
                "unique_distortion_site_count": unique_point_count,
                "embedding_quality": quality,
                "displacement_summary": displacement_summary,
                "verified_sg_num": verified_sg_num,
                "verified_sg_symbol": verified_sg_symbol,
                "predicted_sg": predicted_sg,
                "reference_origin_mode": origin_mode,
                "structure_file": structure_filename,
                "plot_file": plot_filename,
                "signature": signature,
                "residual_match": residual_match,
                "notes": notes,
            }
        )

    timings["candidate_modeling_seconds"] = round(time.perf_counter() - t0, 3)

    coupled_pairs: list[dict[str, object]] = []
    coupled_subgroups: list[dict[str, object]] = []
    possible_coupled_pair_count = 0
    if args.include_coupled_catalog:
        vector_irreps = [row for row in irrep_screen if row["has_vector_channel"]]
        possible_coupled_pair_count = math.comb(len(vector_irreps), 2) if len(vector_irreps) >= 2 else 0
        t0 = time.perf_counter()
        pair_counter = 1
        subgroup_counter = 1
        for irrep_a, irrep_b in combinations(vector_irreps, 2):
            if args.max_coupled_pairs and len(coupled_pairs) >= args.max_coupled_pairs:
                break
            pair_id = f"CP{pair_counter:03d}"
            pair_counter += 1
            subgroups = get_coupled_subgroups(iso, parent_info.sg_num, irrep_a["irrep"], irrep_b["irrep"])
            subgroup_ids = []
            for subgroup in subgroups:
                coupled_id = f"CQ{subgroup_counter:04d}"
                subgroup_counter += 1
                subgroup_ids.append(coupled_id)
                coupled_subgroups.append(
                    {
                        "coupled_id": coupled_id,
                        "pair_id": pair_id,
                        "irrep_a": irrep_a["irrep"],
                        "kpoint_a": irrep_a["kpoint"],
                        "irrep_b": irrep_b["irrep"],
                        "kpoint_b": irrep_b["kpoint"],
                        "reducible_irrep": subgroup["reducible_irrep"],
                        "sg_num": subgroup["sg_num"],
                        "symbol": subgroup["symbol"],
                        "direction": subgroup["direction"],
                        "order_parameter": subgroup["order_parameter"],
                        "basis_text": subgroup["basis_text"],
                        "origin_text": subgroup["origin_text"],
                    }
                )
            coupled_pairs.append(
                {
                    "pair_id": pair_id,
                    "irrep_a": irrep_a["irrep"],
                    "kpoint_a": irrep_a["kpoint"],
                    "irrep_b": irrep_b["irrep"],
                    "kpoint_b": irrep_b["kpoint"],
                    "screen_id_a": irrep_a["screen_id"],
                    "screen_id_b": irrep_b["screen_id"],
                    "subgroup_count": len(subgroups),
                    "subgroup_ids": subgroup_ids,
                }
            )
        timings["coupled_catalog_seconds"] = round(time.perf_counter() - t0, 3)

    unique_subgroups = unique_subgroup_summary(branches, domain_counts)
    signature_rows = []
    for group_index, (key, candidate_ids) in enumerate(sorted(signature_groups.items()), start=1):
        exemplar = next(candidate for candidate in modeled_candidates if candidate.get("signature", {}).get("signature_key") == key)
        signature_rows.append(
            {
                "signature_group_id": f"S{group_index:03d}",
                "candidate_ids": candidate_ids,
                "emergent_peaks": exemplar["signature"]["emergent_peaks"],
                "changed_parent_peaks": exemplar["signature"]["changed_parent_peaks"],
                "suppressed_parent_peaks": exemplar["signature"]["suppressed_parent_peaks"],
            }
        )

    residual_ranking = []
    if residual_peaks:
        residual_ranking = sorted(
            [candidate for candidate in modeled_candidates if candidate.get("residual_match")],
            key=lambda item: item["residual_match"]["score"],
            reverse=True,
        )

    timings["total_seconds"] = round(time.perf_counter() - start_total, 3)

    single_irrep_verified = len([candidate for candidate in modeled_candidates if candidate["status"] == "verified"])
    stats = {
        "iso_call_count": iso.call_count,
        "fixed_kpoint_count": len(fixed_kpoints),
        "parameterized_kpoint_count": len([kpoint for kpoint in kpoints if not kpoint.is_fixed]),
        "catalog_branch_count": len(branches),
        "unique_subgroup_count": len(unique_subgroups),
        "screened_irrep_count": len(irrep_screen),
        "vector_active_irrep_count": len([row for row in irrep_screen if row["has_vector_channel"]]),
        "primary_branch_count": len([branch for branch in branches if branch.is_primary_direction]),
        "modeled_candidate_count": len(modeled_candidates),
        "verified_candidate_count": single_irrep_verified,
        "tentative_candidate_count": len([candidate for candidate in modeled_candidates if candidate["status"] == "tentative"]),
        "embedding_failed_count": len([candidate for candidate in modeled_candidates if candidate["status"] == "embedding_failed"]),
        "signature_group_count": len(signature_rows),
        "coupled_pair_count": len(coupled_pairs),
        "coupled_subgroup_count": len(coupled_subgroups),
        "possible_coupled_pair_count": possible_coupled_pair_count,
    }

    return {
        "settings": {
            "parent_cif": str(parent_path),
            "label": args.label,
            "amplitude": args.amplitude,
            "modeled_direction_prefix": MODELING_DIRECTION_PREFIX,
            "max_modeled": args.max_modeled,
            "include_coupled_catalog": args.include_coupled_catalog,
            "max_coupled_pairs": args.max_coupled_pairs,
            "fixed_kpoints_only": True,
            "parameterized_k_manifolds_documented_only": True,
            "residual_peaks_file": args.residual_peaks,
        },
        "timings": timings,
        "stats": stats,
        "transition_classes": transition_class_rows(
            {
                "include_coupled_catalog": args.include_coupled_catalog,
            },
            stats,
        ),
        "parent": asdict(parent_info),
        "kpoints": [asdict(kpoint) for kpoint in kpoints],
        "branches": [asdict(branch) for branch in branches],
        "irrep_screen": irrep_screen,
        "unique_subgroups": unique_subgroups,
        "modeled_candidates": modeled_candidates,
        "coupled_pairs": coupled_pairs,
        "coupled_subgroups": coupled_subgroups,
        "signature_groups": signature_rows,
        "residual_peaks": residual_peaks,
        "residual_ranking": [
            {
                "candidate_id": candidate["candidate_id"],
                "subgroup_id": candidate["subgroup_id"],
                "score": candidate["residual_match"]["score"],
                "coverage": candidate["residual_match"]["coverage"],
                "precision": candidate["residual_match"]["precision"],
            }
            for candidate in residual_ranking
        ],
    }


def html_escape(text: object) -> str:
    value = str(text)
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def make_html_table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html_escape(header)}</th>" for header in headers)
    body = "".join("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def render_html(report: dict[str, object], run_dir: Path) -> str:
    parent = report["parent"]
    stats = report["stats"]
    summary_rows = [
        ["Parent formula", html_escape(parent["formula"])],
        ["Parent space group", f"{parent['sg_num']} {html_escape(parent['sg_symbol'])}"],
        ["Cataloged single-irrep branches", str(stats["catalog_branch_count"])],
        ["Unique subgroup hypotheses", str(stats["unique_subgroup_count"])],
        ["Screened fixed-k irreps", str(stats["screened_irrep_count"])],
        ["Vector-active irreps", str(stats["vector_active_irrep_count"])],
        ["Primary branches modeled", str(stats["modeled_candidate_count"])],
        ["Symmetry-verified probe candidates", str(stats["verified_candidate_count"])],
        ["Tentative subgroup probes", str(stats["tentative_candidate_count"])],
        ["Embedding failures", str(stats["embedding_failed_count"])],
        ["Coupled irrep pairs cataloged", f"{stats['coupled_pair_count']} / {stats['possible_coupled_pair_count']}"],
        ["Coupled subgroup branches", str(stats["coupled_subgroup_count"])],
        ["Signature groups", str(stats["signature_group_count"])],
        ["ISO calls", str(stats["iso_call_count"])],
    ]

    timing_rows = [[key.replace("_", " "), f"{value:.3f} s"] for key, value in report["timings"].items()]
    subgroup_rows = [
        [
            row["subgroup_id"],
            f"{row['sg_num']} {html_escape(row['symbol'])}",
            str(row["branch_count"]),
            str(row["primary_branch_count"]),
            html_escape(row["basis_text"]),
            html_escape(row["origin_text"]),
            html_escape(row["domain_count"] or "not computed"),
        ]
        for row in report["unique_subgroups"]
    ]
    transition_rows = [
        [html_escape(row["class"]), html_escape(row["status"]), html_escape(row["detail"])]
        for row in report["transition_classes"]
    ]
    irrep_rows = [
        [
            row["screen_id"],
            html_escape(row["kpoint"]),
            html_escape(row["irrep"]),
            str(row["primary_branch_count"]),
            html_escape(row["screening_status"]),
            html_escape(row["exemplar_direction"] or "none"),
            str(row["sample_row_count"]),
        ]
        for row in report["irrep_screen"]
    ]
    coupled_pair_rows = [
        [
            row["pair_id"],
            html_escape(f"{row['kpoint_a']} / {row['irrep_a']}"),
            html_escape(f"{row['kpoint_b']} / {row['irrep_b']}"),
            str(row["subgroup_count"]),
        ]
        for row in report["coupled_pairs"]
    ]
    coupled_subgroup_rows = [
        [
            row["coupled_id"],
            row["pair_id"],
            html_escape(f"{row['sg_num']} {row['symbol']}"),
            html_escape(row["direction"]),
            html_escape(row["order_parameter"]),
        ]
        for row in report["coupled_subgroups"][:60]
    ]

    candidate_rows = []
    detail_blocks = []
    for candidate in report["modeled_candidates"]:
        branch = candidate["branch"]
        signature = candidate.get("signature", {})
        emergent = ", ".join(f"{peak['two_theta']:.2f}" for peak in signature.get("emergent_peaks", [])[:4]) or "none"
        candidate_rows.append(
            [
                candidate["candidate_id"],
                candidate["subgroup_id"],
                html_escape(branch["kpoint"]),
                html_escape(branch["irrep"]),
                html_escape(branch["direction"]),
                html_escape(candidate["status"]),
                html_escape(candidate.get("embedding_quality", "n/a")),
                html_escape(candidate.get("predicted_sg", f"{branch['sg_num']} {branch['symbol']}")),
                emergent,
            ]
        )
        if "signature" not in candidate:
            continue
        residual_html = ""
        if candidate.get("residual_match"):
            residual_html = (
                f"<p><strong>Residual compatibility</strong>: score {candidate['residual_match']['score']:.3f}, "
                f"coverage {candidate['residual_match']['coverage']:.3f}, precision {candidate['residual_match']['precision']:.3f}</p>"
            )
        note_items = "".join(f"<li>{html_escape(note)}</li>" for note in candidate["notes"])
        detail_blocks.append(
            f"""
            <details class="candidate-card" open>
              <summary>
                <span>{html_escape(candidate['candidate_id'])}: {html_escape(branch['kpoint'])} / {html_escape(branch['irrep'])} / {html_escape(branch['direction'])}</span>
                <span class="status-{html_escape(candidate['status'])}">{html_escape(candidate['status'])}</span>
              </summary>
              <div class="candidate-grid">
                <div><strong>Predicted subgroup</strong><br>{html_escape(candidate['predicted_sg'])}</div>
                <div><strong>Verified subgroup</strong><br>{candidate['verified_sg_num']} {html_escape(candidate['verified_sg_symbol'])}</div>
                <div><strong>Basis</strong><br>{html_escape(branch['basis_text'])}</div>
                <div><strong>Origin</strong><br>{html_escape(branch['origin_text'])}</div>
                <div><strong>Reference embedding</strong><br>{html_escape(candidate['reference_origin_mode'])}</div>
                <div><strong>Embedding quality</strong><br>{html_escape(candidate['embedding_quality'])}</div>
                <div><strong>Moved sites</strong><br>{candidate['displacement_summary']['moved_site_count']} ({html_escape(', '.join(candidate['displacement_summary']['moved_species']) or 'none')})</div>
                <div><strong>Match quality</strong><br>{candidate['match']['matched_count']} sites, max distance {candidate['match']['max_distance']:.4f}</div>
                <div><strong>Row collapse</strong><br>{candidate['raw_distortion_row_count']} raw rows to {candidate['unique_distortion_site_count']} unique sites</div>
                <div><strong>Cell size estimate</strong><br>{candidate['reference_build']['actual_site_count']} sites / det {candidate['reference_build']['determinant']:.3f}</div>
              </div>
              {residual_html}
              <p><strong>Emergent peaks</strong>: {html_escape(', '.join(f"{peak['two_theta']:.2f}" for peak in candidate['signature']['emergent_peaks']) or 'none')}</p>
              <p><strong>Suppressed parent peaks</strong>: {html_escape(', '.join(f"{peak['two_theta']:.2f}" for peak in candidate['signature']['suppressed_parent_peaks']) or 'none')}</p>
              <p><strong>Changed parent peaks</strong>: {html_escape(', '.join(f"{peak['two_theta']:.2f}" for peak in candidate['signature']['changed_parent_peaks']) or 'none')}</p>
              <p><strong>Files</strong>: <a href="candidate_structures/{html_escape(candidate['structure_file'])}">{html_escape(candidate['structure_file'])}</a>,
                 <a href="signature_plots/{html_escape(candidate['plot_file'])}">{html_escape(candidate['plot_file'])}</a></p>
              <div class="figure-block">
                <img src="signature_plots/{html_escape(candidate['plot_file'])}" alt="{html_escape(candidate['candidate_id'])}">
              </div>
              <h4>Notes</h4>
              <ul>{note_items or '<li>No additional notes.</li>'}</ul>
            </details>
            """
        )

    signature_rows = [
        [
            row["signature_group_id"],
            html_escape(", ".join(row["candidate_ids"])),
            html_escape(", ".join(f"{peak['two_theta']:.2f}" for peak in row["emergent_peaks"]) or "none"),
            html_escape(", ".join(f"{peak['two_theta']:.2f}" for peak in row["changed_parent_peaks"]) or "none"),
        ]
        for row in report["signature_groups"]
    ]

    residual_rows = [
        [row["candidate_id"], row["subgroup_id"], f"{row['score']:.3f}", f"{row['coverage']:.3f}", f"{row['precision']:.3f}"]
        for row in report["residual_ranking"]
    ]

    methodology = """
    <ol>
      <li>Standardize the parent CIF into a conventional cell.</li>
      <li>Enumerate all fixed special k-points, their irreps, and all single-irrep isotropy-subgroup branches from the local <code>iso</code> executable.</li>
      <li>Group branches into unique subgroup hypotheses using subgroup type, basis, and origin.</li>
      <li>For primary <code>P...</code> directions only, transform <code>DISPLAY DISTORTION</code> points and vectors from the parent-basis display cell into subgroup fractional coordinates, then collapse repeated projected-vector contributions onto unique child-cell sites.</li>
      <li>Construct a canonical subgroup cell, apply a small probe distortion, verify the resulting symmetry, and classify the embedding before using the candidate as evidence.</li>
      <li>Cluster probe candidates by their predicted diffraction signatures and optionally rank them against an external residual-peak list.</li>
    </ol>
    """

    limitations = """
    <ul>
      <li>The catalog is rigorous for fixed-k single-irrep branches. Parameterized k-lines, k-planes, and general q-vectors are documented but not expanded automatically in this first-pass workflow.</li>
      <li>The probe structures are symmetry-guided small-amplitude distortions, not refined physical child structures.</li>
      <li>Coupled-irrep transitions are not searched exhaustively yet; they should be opened only after the single-irrep catalog has been pruned by signatures or data.</li>
      <li>The subgroup cell is reconstructed from basis, origin, and microscopic vectors. This is strong enough for a first-pass signature map, but it is not identical to the full web ISODISTORT child-structure generator.</li>
      <li>Only candidates with both subgroup agreement and exact or good geometric embedding are counted as verified. Tentative and embedding-failed cases stay in the report as auditable dead ends, not as solutions.</li>
    </ul>
    """

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ISOTROPY Discovery Report</title>
  <style>
    :root {{
      --bg: #f5f3ed;
      --panel: #fffdf8;
      --ink: #1f1f1b;
      --muted: #6b665f;
      --accent: #1a6a77;
      --line: #d8d2c8;
      --warn: #935d16;
      --ok: #226a3a;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      background: linear-gradient(180deg, #f9f6f0 0%, var(--bg) 100%);
      line-height: 1.55;
    }}
    main {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 28px 20px 60px;
    }}
    h1, h2, h3, h4 {{
      font-family: "Iowan Old Style", "Palatino Linotype", serif;
      margin: 0 0 0.8rem;
    }}
    h1 {{ font-size: clamp(2rem, 4vw, 3.2rem); }}
    .hero, .panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 1rem 1.1rem;
      box-shadow: 0 12px 28px rgba(45, 33, 24, 0.05);
    }}
    .hero {{
      background: linear-gradient(135deg, rgba(26,106,119,0.10), rgba(147,93,22,0.08));
      margin-bottom: 1.2rem;
    }}
    .grid {{
      display: grid;
      gap: 1rem;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      margin: 1rem 0;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: white;
      border: 1px solid var(--line);
      font-size: 0.94rem;
    }}
    th, td {{
      padding: 0.6rem 0.55rem;
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: top;
    }}
    th {{ background: #f1ebe0; }}
    code {{
      font-family: "IBM Plex Mono", "SFMono-Regular", monospace;
      font-size: 0.92em;
      background: rgba(26,106,119,0.08);
      padding: 0.1em 0.3em;
      border-radius: 4px;
    }}
    .candidate-card {{
      margin: 1rem 0;
      background: rgba(255,255,255,0.9);
      border: 1px solid var(--line);
      border-radius: 14px;
      overflow: hidden;
    }}
    .candidate-card summary {{
      display: flex;
      justify-content: space-between;
      gap: 1rem;
      padding: 0.9rem 1rem;
      cursor: pointer;
      font-weight: 700;
      background: #f7f2e9;
      list-style: none;
    }}
    .candidate-card summary::-webkit-details-marker {{ display: none; }}
    .candidate-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
      gap: 0.9rem;
      padding: 1rem;
    }}
    .candidate-card p, .candidate-card h4, .candidate-card ul {{ padding: 0 1rem; }}
    .figure-block {{ padding: 0 1rem 1rem; }}
    .figure-block img {{
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: white;
    }}
    .status-verified {{ color: var(--ok); }}
    .status-tentative {{ color: var(--warn); }}
    .status-mismatch, .status-no_vector_mode, .status-embedding_failed {{ color: var(--warn); }}
    @media (max-width: 720px) {{
      main {{ padding: 16px 12px 40px; }}
    }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>ISOTROPY Distortion Discovery Report</h1>
      <p>This run treats the parent CIF as a symmetry-search starting point. It catalogs all fixed-k single-irrep subgroup branches locally with <code>iso</code>, then probes the primary branches with small canonical distortions to build a first-pass diffraction signature map.</p>
      {make_html_table(["Metric", "Value"], summary_rows)}
    </section>

    <div class="grid">
      <section class="panel">
        <h2>Search Envelope</h2>
        <p><strong>Parent CIF</strong>: {html_escape(report['settings']['parent_cif'])}</p>
        <p><strong>Probe amplitude</strong>: {report['settings']['amplitude']:.3f}</p>
        <p><strong>Expanded rigorously</strong>: fixed special k-points, single irreps, all subgroup directions</p>
        <p><strong>Modeled expensively</strong>: primary <code>{MODELING_DIRECTION_PREFIX}...</code> directions only</p>
        <p><strong>Coupled-primary catalog</strong>: {'enabled' if report['settings']['include_coupled_catalog'] else 'disabled'}</p>
        <p><strong>Parameterized manifolds</strong>: documented only in this first pass</p>
        <h3>Method</h3>
        {methodology}
      </section>
      <section class="panel">
        <h2>Timing</h2>
        {make_html_table(["Step", "Time"], timing_rows)}
        <h3>Limits</h3>
        {limitations}
      </section>
    </div>

    <section class="panel">
      <h2>Transition-Class Coverage</h2>
      <p>This table is the guardrail against overclaiming. It separates what this run actually searched from what the wider ISOTROPY suite can in principle treat but this local workflow has not yet automated.</p>
      {make_html_table(["Transition class", "Status", "Detail"], transition_rows)}
    </section>

    <section class="panel">
      <h2>Parent Standardization</h2>
      {make_html_table(
          ["Field", "Value"],
          [
              ["Source CIF", html_escape(parent["source_cif"])],
              ["Standardized CIF", html_escape(parent["standardized_cif"])],
              ["Raw SG", f"{parent['raw_sg_num']} {html_escape(parent['raw_sg_symbol'])}"],
              ["Standardized SG", f"{parent['sg_num']} {html_escape(parent['sg_symbol'])}"],
              ["Formula", html_escape(parent["formula"])],
              ["Occupied Wyckoffs", html_escape(', '.join(parent['occupied_wyckoffs']))],
          ],
      )}
    </section>

    <section class="panel">
      <h2>Irrep Structural Screen</h2>
      <p>Each fixed-k irrep was screened for a microscopic vector channel on the occupied Wyckoff sites before coupled-primary intersections were considered. This screen tests symmetry allowance, not dynamical instability or energetic preference.</p>
      {make_html_table(["Screen", "k", "Irrep", "Primary branches", "Screening status", "Example direction", "Sample rows"], irrep_rows)}
    </section>

    <section class="panel">
      <h2>Unique Subgroup Hypotheses</h2>
      <p>These rows collapse multiple branch labels onto subgroup type + basis + origin. This is the right level for reporting the coarse search space before mode-level signatures are compared.</p>
      {make_html_table(["ID", "Subgroup", "Branches", "Primary", "Basis", "Origin", "Domains"], subgroup_rows)}
    </section>

    <section class="panel">
      <h2>Coupled Order Parameters</h2>
      <p>Coupled-primary subgroups are intersections of two uncoupled isotropy subgroups. They matter when no single primary irrep explains an observed child symmetry. If the pair count is capped, the current script expands pairs in the sorted order of the screened irrep table.</p>
      <p><strong>Expanded pairs</strong>: {stats['coupled_pair_count']} of {stats['possible_coupled_pair_count']} screened vector-active irrep pairs.</p>
      {make_html_table(["Pair", "Irrep A", "Irrep B", "Subgroups"], coupled_pair_rows or [["Not run", "", "", ""]])}
      {make_html_table(["Coupled ID", "Pair", "Subgroup", "Direction", "Order parameter"], coupled_subgroup_rows or [["None", "", "", "", ""]])}
    </section>

    <section class="panel">
      <h2>Modeled Primary Branches</h2>
      <p>Only primary branches were carried through to canonical probe structures. The top-level table below is the quick signature index; the expandable cards contain the actual probe details and plots.</p>
      {make_html_table(["Candidate", "Subgroup", "k", "Irrep", "Dir", "Status", "Embedding", "Predicted SG", "Emergent peaks (2θ)"], candidate_rows)}
      {''.join(detail_blocks) or '<p>No primary branches were modeled.</p>'}
    </section>

    <section class="panel">
      <h2>Signature Groups</h2>
      <p>Candidates grouped here have the same coarse diffraction-signature key under the current peak-rounding tolerance.</p>
      {make_html_table(["Group", "Candidates", "Emergent peaks", "Changed parent peaks"], signature_rows or [["None", "", "", ""]])}
    </section>

    <section class="panel">
      <h2>Residual Matching</h2>
      <p>If a residual-peak file is supplied, candidates are ranked by overlap with predicted emergent and changed-peak positions.</p>
      {make_html_table(["Candidate", "Subgroup", "Score", "Coverage", "Precision"], residual_rows or [["No residual file", "", "", "", ""]])}
    </section>
  </main>
</body>
</html>
"""


def render_markdown(report: dict[str, object]) -> str:
    lines = [
        "# ISOTROPY Distortion Discovery Report",
        "",
        "## Summary",
        f"- Parent: `{report['parent']['formula']}` in `{report['parent']['sg_num']} {report['parent']['sg_symbol']}`",
        f"- Cataloged single-irrep branches: `{report['stats']['catalog_branch_count']}`",
        f"- Unique subgroup hypotheses: `{report['stats']['unique_subgroup_count']}`",
        f"- Screened fixed-k irreps: `{report['stats']['screened_irrep_count']}`",
        f"- Vector-active irreps: `{report['stats']['vector_active_irrep_count']}`",
        f"- Modeled primary branches: `{report['stats']['modeled_candidate_count']}`",
        f"- Verified probe candidates: `{report['stats']['verified_candidate_count']}`",
        f"- Tentative probe candidates: `{report['stats']['tentative_candidate_count']}`",
        f"- Embedding failures: `{report['stats']['embedding_failed_count']}`",
        f"- Coupled irrep pairs cataloged: `{report['stats']['coupled_pair_count']} / {report['stats']['possible_coupled_pair_count']}`",
        f"- Coupled subgroup branches: `{report['stats']['coupled_subgroup_count']}`",
        f"- Signature groups: `{report['stats']['signature_group_count']}`",
        f"- ISO calls: `{report['stats']['iso_call_count']}`",
        "",
        "## Scope",
        "- Fixed special k-points are expanded rigorously.",
        "- Parameterized k-lines/planes/general q vectors are documented only in this first pass.",
        "- Primary `P...` branches are the only ones carried through to canonical probe distortions.",
        f"- Coupled-primary catalog: `{'enabled' if report['settings']['include_coupled_catalog'] else 'disabled'}`.",
        "- `DISPLAY DISTORTION` rows are transformed from the parent-basis display cell into subgroup fractional coordinates before site matching.",
        "- Repeated projected-vector contributions on the same child site are summed before the probe distortion is applied.",
        "- Probe distortions are signature probes, not refined child structures.",
        "",
        "## Transition-Class Coverage",
    ]
    for row in report["transition_classes"]:
        lines.append(f"- `{row['class']}`: `{row['status']}`. {row['detail']}")
    lines.extend(
        [
            "",
            "## Timings",
        ]
    )
    lines.extend(f"- `{key}`: {value:.3f} s" for key, value in report["timings"].items())
    lines.extend(
        [
            "",
            "## Irrep Structural Screen",
        ]
    )
    lines.append("- `vector_active` means symmetry allows a microscopic displacement channel on the occupied parent sites. It does not by itself imply a soft mode or an energetically preferred transition.")
    for row in report["irrep_screen"]:
        lines.append(
            f"- `{row['screen_id']}` `{row['kpoint']} / {row['irrep']}` -> `{row['screening_status']}`; primary branches `{row['primary_branch_count']}`, example direction `{row['exemplar_direction'] or 'none'}`"
        )
    lines.extend(
        [
            "",
            "## Unique Subgroups",
        ]
    )
    for row in report["unique_subgroups"]:
        lines.append(
            f"- `{row['subgroup_id']}`: `{row['sg_num']} {row['symbol']}` from {row['branch_count']} branches; basis `{row['basis_text']}`, origin `{row['origin_text']}`"
        )
    lines.extend(
        [
            "",
            "## Coupled Order Parameters",
        ]
    )
    for row in report["coupled_pairs"]:
        lines.append(
            f"- `{row['pair_id']}` `{row['kpoint_a']} / {row['irrep_a']}` + `{row['kpoint_b']} / {row['irrep_b']}` -> `{row['subgroup_count']}` coupled subgroup branches"
        )
    if report["settings"]["include_coupled_catalog"]:
        lines.append(
            f"- Expanded `{report['stats']['coupled_pair_count']}` of `{report['stats']['possible_coupled_pair_count']}` screened vector-active irrep pairs under the current cap."
        )
    lines.extend(
        [
            "",
            "## Modeled Candidates",
        ]
    )
    for candidate in report["modeled_candidates"]:
        branch = candidate["branch"]
        line = (
            f"- `{candidate['candidate_id']}` `{branch['kpoint']} / {branch['irrep']} / {branch['direction']}` "
            f"-> `{candidate['status']}` with embedding `{candidate.get('embedding_quality', 'n/a')}`"
        )
        if candidate.get("signature"):
            emergent = ", ".join(f"{peak['two_theta']:.2f}" for peak in candidate["signature"]["emergent_peaks"]) or "none"
            line += f"; emergent peaks `{emergent}`"
        lines.append(line)
    if report["residual_ranking"]:
        lines.extend(
            [
                "",
                "## Residual Ranking",
            ]
        )
        for row in report["residual_ranking"]:
            lines.append(f"- `{row['candidate_id']}` score `{row['score']:.3f}` coverage `{row['coverage']:.3f}` precision `{row['precision']:.3f}`")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Catalog symmetry-breaking branches and build canonical diffraction signatures.")
    parser.add_argument("parent_cif", help="Parent CIF file.")
    parser.add_argument("--label", help="Output run label. Defaults to the parent stem.")
    parser.add_argument("--amplitude", type=float, default=PROBE_AMPLITUDE, help="Canonical probe amplitude in fractional mode units.")
    parser.add_argument("--max-modeled", type=int, default=60, help="Maximum number of primary branches to model.")
    parser.add_argument("--include-coupled-catalog", action="store_true", help="Also catalog coupled-primary fixed-k subgroup intersections for screened vector-active irreps.")
    parser.add_argument("--max-coupled-pairs", type=int, default=20, help="Maximum number of screened irrep pairs to expand in the coupled-primary catalog.")
    parser.add_argument("--residual-peaks", help="Optional CSV/JSON file of residual peak positions for ranking candidates.")
    args = parser.parse_args()

    if not args.label:
        args.label = Path(args.parent_cif).stem.replace(".", "_")

    report = run_pipeline(args)
    run_dir = DISCOVERY_RUNS_DIR / args.label

    (run_dir / "discovery_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    (run_dir / "discovery_report.md").write_text(render_markdown(report), encoding="utf-8")
    (run_dir / "discovery_report.html").write_text(render_html(report, run_dir), encoding="utf-8")


if __name__ == "__main__":
    main()

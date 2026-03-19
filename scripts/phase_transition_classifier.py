from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from _isotropy_paths import configure_environment, resolve_input_path
from isodistort_distortion_parser import parse_distortion_file


configure_environment()


@dataclass
class TransitionClassification:
    family: str
    confidence: str
    rationale: list[str]
    parent_path: str | None
    child_path: str | None
    distortion_path: str | None
    parent_space_group: str | None
    child_space_group: str | None
    recommended_workflow: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _structure_and_sg(path: Path) -> tuple[object, str]:
    structure = CifParser(str(path)).parse_structures(primitive=False)[0]
    analyzer = SpacegroupAnalyzer(structure, symprec=1e-2)
    return structure, f"{analyzer.get_space_group_number()} {analyzer.get_space_group_symbol()}"


def classify_transition(
    *,
    parent_cif: str | None = None,
    child_cif: str | None = None,
    distortion_file: str | None = None,
) -> TransitionClassification:
    parent_path = resolve_input_path(parent_cif) if parent_cif else None
    child_path = resolve_input_path(child_cif) if child_cif else None
    distortion_path = resolve_input_path(distortion_file) if distortion_file else None

    if distortion_path is not None:
        distortion_data = parse_distortion_file(distortion_path)
        rationale = [
            "An ISODISTORT distortion file was supplied, so the most direct workflow is explanation/decomposition rather than blind parent-side discovery."
        ]
        workflow = "generate_tutorial_reports_or_distortion_explainer"
        family = "decomposition_from_distortion_file"
        if parent_path is not None:
            rationale.append("A parent CIF is also available, so the explanation can anchor the distortion in the parent symmetry setting.")
        if child_path is not None:
            rationale.append("A child CIF is also available, so the workflow can check whether parent and known child settings align with the saved decomposition.")
        if child_path is not None and "magnetic" in distortion_data.classification:
            rationale.append("The distortion file carries magnetic modes and a known child structure was supplied, so the dedicated magnetic parent/child workflow is the best fit.")
            workflow = "analyze_magnetic_child_structure"
            family = "known_magnetic_parent_child_decomposition"
        return TransitionClassification(
            family=family,
            confidence="high",
            rationale=rationale,
            parent_path=str(parent_path) if parent_path else None,
            child_path=str(child_path) if child_path else None,
            distortion_path=str(distortion_path),
            parent_space_group=None,
            child_space_group=None,
            recommended_workflow=workflow,
        )

    if parent_path is not None and child_path is None:
        return TransitionClassification(
            family="symmetry_connected_parent_only",
            confidence="high",
            rationale=["Only a parent structure was supplied, so the correct starting point is parent-side discovery over irreps and isotropy subgroups."],
            parent_path=str(parent_path),
            child_path=None,
            distortion_path=None,
            parent_space_group=None,
            child_space_group=None,
            recommended_workflow="discover_distortion_signatures",
        )

    if parent_path is None and child_path is not None:
        return TransitionClassification(
            family="child_only_unsupported",
            confidence="high",
            rationale=["A child structure without a parent is insufficient for the intended ISOTROPY workflows in this repository."],
            parent_path=None,
            child_path=str(child_path),
            distortion_path=None,
            parent_space_group=None,
            child_space_group=None,
            recommended_workflow="supply_parent_structure",
        )

    if parent_path is None and child_path is None:
        return TransitionClassification(
            family="insufficient_input",
            confidence="high",
            rationale=["No parent CIF, child CIF, or distortion file was supplied."],
            parent_path=None,
            child_path=None,
            distortion_path=None,
            parent_space_group=None,
            child_space_group=None,
            recommended_workflow="supply_parent_structure_or_distortion_file",
        )

    parent_structure, parent_sg = _structure_and_sg(parent_path)
    child_structure, child_sg = _structure_and_sg(child_path)
    matcher = StructureMatcher(primitive_cell=False, scale=True, attempt_supercell=True)
    same_framework = matcher.fit(parent_structure, child_structure)

    rationale = [
        f"Parent structure analyzed as {parent_sg}.",
        f"Child structure analyzed as {child_sg}.",
    ]

    if same_framework:
        rationale.append("The parent and child are structurally matchable under a tolerant structure matcher, suggesting a symmetry-connected workflow is plausible.")
        family = "symmetry_connected_parent_child"
        confidence = "medium"
        workflow = "discover_then_compare_or_mode_decompose"
    else:
        rationale.append("The parent and child are not directly matchable under a tolerant structure matcher, so a reconstructive/common-subgroup workflow is the safer default.")
        family = "reconstructive_parent_child"
        confidence = "medium"
        workflow = "reconstructive_common_subgroups"

    return TransitionClassification(
        family=family,
        confidence=confidence,
        rationale=rationale,
        parent_path=str(parent_path),
        child_path=str(child_path),
        distortion_path=None,
        parent_space_group=parent_sg,
        child_space_group=child_sg,
        recommended_workflow=workflow,
    )

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from _isotropy_paths import configure_environment, resolve_input_path
from transition_metadata import inspect_cif_metadata, inspect_distortion_metadata


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
    detected_metadata: dict[str, dict[str, object]]

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

    detected_metadata: dict[str, dict[str, object]] = {}
    if parent_path is not None:
        detected_metadata["parent"] = inspect_cif_metadata(parent_path).as_dict()
    if child_path is not None:
        detected_metadata["child"] = inspect_cif_metadata(child_path).as_dict()
    if distortion_path is not None:
        detected_metadata["distortion"] = inspect_distortion_metadata(distortion_path).as_dict()

    magnetic_present = any(item.get("has_magnetic_metadata") for item in detected_metadata.values())
    superspace_present = any(item.get("has_superspace_metadata") for item in detected_metadata.values())
    coupled_structural_present = any(
        item.get("kind") == "distortion" and item.get("has_structural_distortion_metadata") for item in detected_metadata.values()
    )

    if distortion_path is not None:
        rationale = [
            "An ISODISTORT distortion file was supplied, so routing should respect the distortion metadata rather than defaulting to generic structural reporting."
        ]
        if superspace_present and magnetic_present:
            rationale.append("The supplied distortion metadata indicate an incommensurate magnetic superspace case.")
            family = "incommensurate_magnetic_superspace_analysis"
            workflow = "incommensurate_magnetic_superspace_analysis"
        elif magnetic_present and coupled_structural_present:
            rationale.append("The supplied distortion metadata include both magnetic and structural channels, so coupled analysis is the right family.")
            family = "mixed_magnetic_structural_coupled_analysis"
            workflow = "mixed_magnetic_structural_coupled_analysis"
        elif magnetic_present:
            rationale.append("The supplied distortion metadata indicate a magnetic decomposition problem.")
            family = "magnetic_parent_child_decomposition"
            workflow = "magnetic_parent_child_decomposition"
        else:
            if parent_path is not None:
                rationale.append("A parent CIF is also available, so the explanation can anchor the distortion in the parent symmetry setting.")
            family = "decomposition_from_distortion_file"
            workflow = "generate_tutorial_reports_or_distortion_explainer"
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
            detected_metadata=detected_metadata,
        )

    if parent_path is not None and child_path is None:
        rationale = ["Only a parent structure was supplied, so the router should start from parent-side discovery metadata."]
        if superspace_present and magnetic_present:
            rationale.append("The supplied parent CIF already carries magnetic superspace tags, so it belongs to a dedicated superspace workflow.")
            family = "incommensurate_magnetic_superspace_analysis"
            workflow = "incommensurate_magnetic_superspace_analysis"
        elif magnetic_present:
            rationale.append("The supplied parent CIF contains magnetic metadata, so it should route to commensurate magnetic parent-only discovery instead of structural discovery.")
            family = "commensurate_magnetic_parent_only_discovery"
            workflow = "commensurate_magnetic_parent_only_discovery"
        else:
            rationale.append("No magnetic or superspace metadata were detected, so structural parent-side discovery remains the default.")
            family = "symmetry_connected_parent_only"
            workflow = "discover_distortion_signatures"
        return TransitionClassification(
            family=family,
            confidence="high",
            rationale=rationale,
            parent_path=str(parent_path),
            child_path=None,
            distortion_path=None,
            parent_space_group=None,
            child_space_group=None,
            recommended_workflow=workflow,
            detected_metadata=detected_metadata,
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
            detected_metadata=detected_metadata,
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
            detected_metadata=detected_metadata,
        )

    parent_structure, parent_sg = _structure_and_sg(parent_path)
    child_structure, child_sg = _structure_and_sg(child_path)
    matcher = StructureMatcher(primitive_cell=False, scale=True, attempt_supercell=True)
    same_framework = matcher.fit(parent_structure, child_structure)

    rationale = [
        f"Parent structure analyzed as {parent_sg}.",
        f"Child structure analyzed as {child_sg}.",
    ]

    if superspace_present and magnetic_present:
        rationale.append("Magnetic superspace metadata were detected in the supplied CIF inputs, so the dedicated superspace workflow takes precedence over structural matching.")
        family = "incommensurate_magnetic_superspace_analysis"
        confidence = "high"
        workflow = "incommensurate_magnetic_superspace_analysis"
    elif magnetic_present and coupled_structural_present:
        rationale.append("Both magnetic and structural metadata were detected across the supplied inputs, so a mixed coupled-analysis workflow is preferred.")
        family = "mixed_magnetic_structural_coupled_analysis"
        confidence = "high"
        workflow = "mixed_magnetic_structural_coupled_analysis"
    elif magnetic_present:
        rationale.append("Magnetic metadata were detected across the supplied inputs, so a dedicated magnetic parent+child decomposition workflow is preferred.")
        family = "magnetic_parent_child_decomposition"
        confidence = "high"
        workflow = "magnetic_parent_child_decomposition"
    elif same_framework:
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
        detected_metadata=detected_metadata,
    )

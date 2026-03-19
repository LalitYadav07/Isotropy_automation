from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import re

from isodistort_distortion_parser import parse_distortion_file


MAGNETIC_CIF_MARKERS = (
    "_space_group_magn",
    "_atom_site_moment",
    "_iso_magneticmode",
    "magnetic space group",
)

SUPERSPACE_CIF_MARKERS = (
    "_cell_modulation_dimension",
    "_cell_wave_vector",
    "_space_group_ssg",
    "_space_group_magn.ssg",
    "_space_group_magn_ssg",
    "superspace group",
)

STRUCTURAL_DISTORTION_CLASSES = {"displacive", "occupational", "rotational", "strain"}


@dataclass
class FileMetadata:
    path: str
    kind: str
    has_magnetic_metadata: bool
    has_superspace_metadata: bool
    has_structural_distortion_metadata: bool
    descriptors: list[str]
    evidence: list[str]

    def as_dict(self) -> dict[str, object]:
        return asdict(self)



def _append_flag(condition: bool, descriptor: str, evidence: str, descriptors: list[str], evidence_items: list[str]) -> None:
    if condition:
        descriptors.append(descriptor)
        evidence_items.append(evidence)



def inspect_cif_metadata(path: Path) -> FileMetadata:
    text = path.read_text(encoding="utf-8", errors="ignore")
    lowered = text.lower()
    descriptors: list[str] = []
    evidence: list[str] = []

    has_magnetic = any(marker in lowered for marker in MAGNETIC_CIF_MARKERS)
    has_superspace = any(marker in lowered for marker in SUPERSPACE_CIF_MARKERS)
    has_structural = bool(re.search(r"_atom_site_(fract|occupancy)|_atom_site_aniso", lowered))

    _append_flag(has_magnetic, "magnetic", "CIF contains magnetic-space-group or magnetic-moment tags.", descriptors, evidence)
    _append_flag(has_superspace, "superspace", "CIF contains modulation-wave or superspace-group tags.", descriptors, evidence)
    _append_flag(has_structural, "structural", "CIF contains conventional atom-site structural data.", descriptors, evidence)

    return FileMetadata(
        path=str(path),
        kind="cif",
        has_magnetic_metadata=has_magnetic,
        has_superspace_metadata=has_superspace,
        has_structural_distortion_metadata=has_structural,
        descriptors=descriptors,
        evidence=evidence,
    )



def inspect_distortion_metadata(path: Path) -> FileMetadata:
    report = parse_distortion_file(path)
    descriptors: list[str] = []
    evidence: list[str] = []

    classes = set(report.classification)
    has_magnetic = "magnetic" in classes
    has_superspace = any(
        [
            report.modulation_count > 0,
            "incommensurate/modulated" in classes,
            "superspace" in report.order_parameter.lower(),
            (report.subgroup_string or "").lower().find("superspace") >= 0,
        ]
    )
    has_structural = bool(STRUCTURAL_DISTORTION_CLASSES & classes)

    _append_flag(has_magnetic, "magnetic", "Distortion file enables magnetic channels or stores magnetic irreps.", descriptors, evidence)
    _append_flag(has_superspace, "superspace", "Distortion file records incommensurate modulation or superspace metadata.", descriptors, evidence)
    _append_flag(has_structural, "structural", "Distortion file enables structural channels such as displacive, occupational, rotational, or strain modes.", descriptors, evidence)

    return FileMetadata(
        path=str(path),
        kind="distortion",
        has_magnetic_metadata=has_magnetic,
        has_superspace_metadata=has_superspace,
        has_structural_distortion_metadata=has_structural,
        descriptors=descriptors,
        evidence=evidence,
    )

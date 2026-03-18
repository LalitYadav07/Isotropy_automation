from __future__ import annotations

import json
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

try:
    import matplotlib.pyplot as plt
    import numpy as np
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.io.cif import CifParser, CifWriter
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
except ImportError as exc:
    raise SystemExit(
        "Missing dependency. Run this script with the project virtualenv, for example:\n"
        "  ./.venv/bin/python isotropy_project/scripts/simulate_diffraction.py"
    ) from exc

try:
    from _isotropy_paths import (
        PLOTS_DIR,
        REPORT_PATH,
        SCRATCH_DIR,
        SUBGROUPS_DIR,
        configure_environment,
        resolve_input_path,
    )
except ModuleNotFoundError:
    from isotropy_project.scripts._isotropy_paths import (
        PLOTS_DIR,
        REPORT_PATH,
        SCRATCH_DIR,
        SUBGROUPS_DIR,
        configure_environment,
        resolve_input_path,
    )


REPORT_HTML_PATH = REPORT_PATH.with_suffix(".html")
REPORT_JSON_PATH = REPORT_PATH.with_suffix(".json")
XRD_WAVELENGTH = "CuKa"
REPORT_2THETA_MIN = 10
REPORT_2THETA_MAX = 90
SYMPREC = 1e-2
DISTORTION_AMPLITUDE = 0.01
SCREEN_WIDTH = 240
PAGE_LENGTH = 1000
ISO_TIMEOUT_SECONDS = 3
MODELED_KPOINTS = {"GM", "L", "X"}


@dataclass
class SiteSummary:
    index: int
    species: str
    wyckoff: str
    frac_coords: list[float]


@dataclass
class ParentInfo:
    cif_path: str
    formula: str
    sg_num: int
    sg_symbol: str
    lattice_parameters: list[float]
    site_summaries: list[SiteSummary]
    wyckoff_letters: list[str]


@dataclass
class KPointInfo:
    label: str
    coordinates: str
    is_fixed: bool


@dataclass
class SubgroupInfo:
    sg_num: int
    symbol: str
    direction: str
    order_parameter: str
    basis_text: str
    basis_matrix: list[list[float]]
    origin_text: str


@dataclass
class DistortionRow:
    wyckoff: str
    point: list[float]
    vector: list[float]
    site_index: int
    species: str


@dataclass
class ModeResult:
    kpoint: str
    irrep: str
    subgroup: SubgroupInfo
    domain_count: int
    displacement_rows: list[DistortionRow]
    displaced_site_count: int
    displaced_species: list[str]
    max_displacement_norm: float
    rms_displacement_norm: float
    verification_status: str
    verified_sg_num: int | None
    verified_sg_symbol: str | None
    xrd_plot_file: str | None
    structure_file: str | None
    xrd_summary: dict[str, list[dict[str, float]]] | None
    notes: list[str]


configure_environment()


def run_iso(commands: str) -> str:
    full_cmd = (
        f"PAGE {PAGE_LENGTH}\n"
        f"SCREEN {SCREEN_WIDTH}\n"
        f"{commands.strip()}\n"
        "QUIT\n"
    )
    try:
        result = subprocess.run(
            ["iso"],
            input=full_cmd,
            capture_output=True,
            text=True,
            cwd=SCRATCH_DIR,
            check=False,
            timeout=ISO_TIMEOUT_SECONDS,
        )
        return result.stdout
    except subprocess.TimeoutExpired as exc:
        return exc.stdout or ""


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


def structure_signature(structure) -> str:
    analyzer = SpacegroupAnalyzer(structure, symprec=SYMPREC)
    symm = analyzer.get_symmetrized_structure()
    species_positions = []
    for site in symm:
        coords = ",".join(f"{value:.5f}" for value in np.mod(site.frac_coords, 1.0))
        species_positions.append(f"{site.species_string}:{coords}")
    species_positions.sort()
    return f"{analyzer.get_space_group_number()}|{analyzer.get_space_group_symbol()}|{'|'.join(species_positions)}"


def get_parent_info(parent_cif: Path):
    structure = CifParser(parent_cif).parse_structures()[0]
    analyzer = SpacegroupAnalyzer(structure, symprec=1e-3)
    dataset = analyzer.get_symmetry_dataset()
    site_summaries: list[SiteSummary] = []
    wyckoff_letters: list[str] = []
    for idx, (site, wyckoff_letter) in enumerate(zip(structure, dataset.wyckoffs, strict=False)):
        wyckoff_letter = wyckoff_letter.lower()
        wyckoff_letters.append(wyckoff_letter)
        site_summaries.append(
            SiteSummary(
                index=idx,
                species=site.species_string,
                wyckoff=wyckoff_letter,
                frac_coords=[float(value) for value in np.mod(site.frac_coords, 1.0)],
            )
        )

    parent_info = ParentInfo(
        cif_path=str(parent_cif),
        formula=structure.composition.reduced_formula,
        sg_num=analyzer.get_space_group_number(),
        sg_symbol=analyzer.get_space_group_symbol(),
        lattice_parameters=[float(value) for value in structure.lattice.abc + structure.lattice.angles],
        site_summaries=site_summaries,
        wyckoff_letters=sorted(set(wyckoff_letters)),
    )
    return parent_info, structure


def get_kpoints(parent_sg: int) -> list[KPointInfo]:
    out = run_iso(
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
        is_fixed = re.search(r"[a-z]", coordinates) is None
        kpoints.append(KPointInfo(match.group(1), coordinates, is_fixed))
    return kpoints


def get_irreps(parent_sg: int, kpoint: str) -> list[str]:
    out = run_iso(
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


def get_subgroups(parent_sg: int, kpoint: str, irrep: str) -> list[SubgroupInfo]:
    out = run_iso(
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
    subgroups: list[SubgroupInfo] = []
    pattern = re.compile(
        r"^(\d+)\s+(\S+)\s+(\S+)\s+(\(.+?\))\s+(\([^)]*\),\([^)]*\),\([^)]*\))\s+(\([^)]*\))$"
    )
    for line in out.splitlines():
        match = pattern.match(line.strip())
        if not match:
            continue
        basis_text = match.group(5)
        subgroups.append(
            SubgroupInfo(
                sg_num=int(match.group(1)),
                symbol=match.group(2),
                direction=match.group(3),
                order_parameter=match.group(4),
                basis_text=basis_text,
                basis_matrix=parse_basis_vectors(basis_text),
                origin_text=match.group(6),
            )
        )
    return subgroups


def get_domain_count(parent_sg: int, kpoint: str, irrep: str, direction: str) -> int:
    out = run_iso(
        f"""
        VALUE PARENT {parent_sg}
        VALUE KPOINT {kpoint}
        VALUE IRREP {irrep}
        VALUE DIRECTION {direction}
        SHOW DOMAINS
        DISPLAY ISOTROPY
        """
    )
    domains = []
    for line in out.splitlines():
        stripped = line.strip()
        if stripped.isdigit():
            domains.append(int(stripped))
    return len(domains)


def get_vector_distortion_rows(
    parent_sg: int,
    kpoint: str,
    irrep: str,
    direction: str,
    wyckoff_letters: list[str],
    domain: int | None = 1,
) -> list[tuple[str, list[float], list[float]]]:
    wyckoff_str = " ".join(letter.upper() for letter in wyckoff_letters)
    domain_cmd = f"VALUE DOMAIN {domain}\n" if domain else ""
    out = run_iso(
        f"""
        VALUE PARENT {parent_sg}
        VALUE KPOINT {kpoint}
        VALUE IRREP {irrep}
        VALUE DIRECTION {direction}
        {domain_cmd}VALUE CELL 1,0,0 0,1,0 0,0,1
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


def site_distance(frac_a: np.ndarray, frac_b: np.ndarray) -> float:
    delta = np.mod(frac_a - frac_b + 0.5, 1.0) - 0.5
    return float(np.linalg.norm(delta))


def build_distorted_structure(parent_structure, parent_info: ParentInfo, rows, amplitude: float):
    distorted = parent_structure.copy()
    applied_rows: list[DistortionRow] = []
    used_indices: set[int] = set()
    site_by_index = {site.index: site for site in parent_info.site_summaries}

    for wyckoff, point, vector in rows:
        point_array = np.mod(np.array(point, dtype=float), 1.0)
        best_match: tuple[float, int] | None = None
        for site in parent_info.site_summaries:
            if site.wyckoff != wyckoff:
                continue
            distance = site_distance(np.array(site.frac_coords), point_array)
            if best_match is None or distance < best_match[0]:
                best_match = (distance, site.index)
        if best_match is None or best_match[0] > 1e-4:
            raise ValueError(f"Could not match distortion point {point} to Wyckoff {wyckoff}.")
        site_index = best_match[1]
        if site_index in used_indices:
            raise ValueError(f"Duplicate displacement mapping for site index {site_index}.")
        used_indices.add(site_index)
        vector_array = np.array(vector, dtype=float)
        new_coords = np.mod(distorted[site_index].frac_coords + amplitude * vector_array, 1.0)
        distorted.replace(site_index, distorted[site_index].species, new_coords, coords_are_cartesian=False)
        applied_rows.append(
            DistortionRow(
                wyckoff=wyckoff,
                point=[float(value) for value in point_array],
                vector=[float(value) for value in vector_array],
                site_index=site_index,
                species=site_by_index[site_index].species,
            )
        )
    return distorted, applied_rows


def summarize_displacements(rows: list[DistortionRow]) -> tuple[int, list[str], float, float]:
    if not rows:
        return 0, [], 0.0, 0.0
    norms = [float(np.linalg.norm(row.vector)) for row in rows]
    species = sorted({row.species for row in rows})
    return len(rows), species, max(norms), float(np.sqrt(np.mean(np.square(norms))))


def verify_structure(structure) -> tuple[int, str]:
    analyzer = SpacegroupAnalyzer(structure, symprec=SYMPREC)
    return analyzer.get_space_group_number(), analyzer.get_space_group_symbol()


def strongest_peaks(pattern, count: int = 6) -> list[dict[str, float]]:
    peaks = list(zip(pattern.x, pattern.y, strict=False))
    peaks.sort(key=lambda item: item[1], reverse=True)
    selected = []
    for two_theta, intensity in peaks:
        if REPORT_2THETA_MIN <= two_theta <= REPORT_2THETA_MAX:
            selected.append({"two_theta": float(two_theta), "intensity": float(intensity)})
        if len(selected) == count:
            break
    return selected


def new_peaks(parent_pattern, distorted_pattern, tolerance: float = 0.15, threshold: float = 5.0):
    parent_peaks = strongest_peaks(parent_pattern, count=20)
    distorted_peaks = strongest_peaks(distorted_pattern, count=20)
    new_peak_list = []
    for peak in distorted_peaks:
        if peak["intensity"] < threshold:
            continue
        if not any(abs(peak["two_theta"] - ref["two_theta"]) <= tolerance for ref in parent_peaks):
            new_peak_list.append(peak)
    return new_peak_list[:6]


def make_xrd_plot(parent_structure, distorted_structure, label: str) -> tuple[str, dict[str, list[dict[str, float]]]]:
    calculator = XRDCalculator(wavelength=XRD_WAVELENGTH)
    parent_pattern = calculator.get_pattern(parent_structure)
    distorted_pattern = calculator.get_pattern(distorted_structure)

    plt.figure(figsize=(10, 6))
    plt.vlines(parent_pattern.x, 0, parent_pattern.y, color="#4c72b0", alpha=0.7, label="Parent")
    plt.vlines(distorted_pattern.x, 0, distorted_pattern.y, color="#dd8452", alpha=0.7, label="Distorted")
    plt.title(label)
    plt.xlabel("2θ (degrees)")
    plt.ylabel("Relative intensity")
    plt.xlim(REPORT_2THETA_MIN, REPORT_2THETA_MAX)
    plt.grid(True, alpha=0.25)
    plt.legend()

    safe_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
    filename = f"{safe_label}.png"
    output_path = PLOTS_DIR / filename
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()

    summary = {
        "parent_strongest": strongest_peaks(parent_pattern),
        "distorted_strongest": strongest_peaks(distorted_pattern),
        "new_peaks": new_peaks(parent_pattern, distorted_pattern),
    }
    return filename, summary


def write_structure_file(structure, name: str) -> str:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    filename = f"{safe_name}.cif"
    path = SUBGROUPS_DIR / filename
    CifWriter(structure).write_file(path)
    return filename


def html_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def make_html_table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html_escape(header)}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def generate_reports(report_data: dict) -> None:
    with open(REPORT_JSON_PATH, "w", encoding="utf-8") as handle:
        json.dump(report_data, handle, indent=2)

    parent = report_data["parent"]
    fixed_rows = [
        [item["label"], item["coordinates"]]
        for item in report_data["kpoints"]
        if item["is_fixed"]
    ]
    parameterized_rows = [
        [item["label"], item["coordinates"]]
        for item in report_data["kpoints"]
        if not item["is_fixed"]
    ]
    site_rows = [
        [
            str(site["index"]),
            site["species"],
            site["wyckoff"],
            ", ".join(f"{value:.3f}" for value in site["frac_coords"]),
        ]
        for site in report_data["parent"]["site_summaries"]
    ]
    catalog_rows = []
    for entry in report_data["catalog"]:
        catalog_rows.append(
            [
                entry["kpoint"],
                entry["irrep"],
                entry["direction"],
                f"{entry['sg_num']} {entry['symbol']}",
                html_escape(entry["order_parameter"]),
                html_escape(entry["basis_text"]),
                html_escape(entry["origin_text"]),
                str(entry["domain_count"]),
            ]
        )

    modeled_sections = []
    for mode in report_data["mode_results"]:
        status_class = {
            "verified": "status-ok",
            "mismatch": "status-warn",
            "no_vector_mode": "status-neutral",
            "construction_failed": "status-warn",
        }.get(mode["verification_status"], "status-neutral")
        note_items = "".join(f"<li>{html_escape(note)}</li>" for note in mode["notes"])
        xrd_block = ""
        if mode["xrd_plot_file"]:
            strongest_rows = [
                [f"{peak['two_theta']:.2f}", f"{peak['intensity']:.1f}"]
                for peak in mode["xrd_summary"]["distorted_strongest"]
            ]
            new_peak_rows = [
                [f"{peak['two_theta']:.2f}", f"{peak['intensity']:.1f}"]
                for peak in mode["xrd_summary"]["new_peaks"]
            ]
            xrd_block = f"""
            <div class="figure-block">
              <img src="diffraction_plots/{html_escape(mode['xrd_plot_file'])}" alt="{html_escape(mode['kpoint'])} {html_escape(mode['irrep'])}">
            </div>
            <div class="two-col">
              <div>
                <h4>Strongest Distorted Peaks</h4>
                {make_html_table(["2θ", "Intensity"], strongest_rows)}
              </div>
              <div>
                <h4>Peaks Not Present In Parent Top Reflections</h4>
                {make_html_table(["2θ", "Intensity"], new_peak_rows or [["None", ""]])}
              </div>
            </div>
            """

        modeled_sections.append(
            f"""
            <details class="mode-card" open>
              <summary>
                <span>{html_escape(mode['kpoint'])} / {html_escape(mode['irrep'])} / {html_escape(mode['subgroup']['direction'])}</span>
                <span class="{status_class}">{html_escape(mode['verification_status'])}</span>
              </summary>
              <div class="mode-grid">
                <div><strong>Predicted subgroup</strong><br>{mode['subgroup']['sg_num']} {html_escape(mode['subgroup']['symbol'])}</div>
                <div><strong>Verified subgroup</strong><br>{mode['verified_sg_num'] or 'N/A'} {html_escape(mode['verified_sg_symbol'] or '')}</div>
                <div><strong>Order-parameter direction</strong><br>{html_escape(mode['subgroup']['order_parameter'])}</div>
                <div><strong>Domain count</strong><br>{mode['domain_count']}</div>
                <div><strong>Basis</strong><br>{html_escape(mode['subgroup']['basis_text'])}</div>
                <div><strong>Origin</strong><br>{html_escape(mode['subgroup']['origin_text'])}</div>
                <div><strong>Displaced sites</strong><br>{mode['displaced_site_count']} ({', '.join(mode['displaced_species']) or 'none'})</div>
                <div><strong>Max / RMS mode norm</strong><br>{mode['max_displacement_norm']:.3f} / {mode['rms_displacement_norm']:.3f}</div>
                <div><strong>Structure file</strong><br>{html_escape(mode['structure_file'] or 'not written')}</div>
              </div>
              <h4>Notes</h4>
              <ul>{note_items}</ul>
              {xrd_block}
            </details>
            """
        )

    summary_rows = [
        ["Parent formula", html_escape(parent["formula"])],
        ["Parent space group", f"{parent['sg_num']} {html_escape(parent['sg_symbol'])}"],
        [
            "Lattice parameters",
            ", ".join(f"{value:.4f}" for value in parent["lattice_parameters"]),
        ],
        ["Fixed k-points analyzed", str(report_data["summary"]["fixed_kpoint_count"])],
        ["Parameterized k-manifolds documented only", str(report_data["summary"]["parameterized_kpoint_count"])],
        ["Subgroup directions cataloged", str(report_data["summary"]["catalog_count"])],
        ["Microscopic vector modes reconstructed", str(report_data["summary"]["mode_count"])],
        ["Symmetry-verified distorted structures", str(report_data["summary"]["verified_mode_count"])],
    ]

    methodology_notes = """
    <ol>
      <li>Identify the parent structure from the input CIF using <code>pymatgen</code>/<code>spglib</code>.</li>
      <li>Query ISOTROPY for fixed k-points, irreps, isotropy subgroups, basis vectors, origins, and domain multiplicities.</li>
      <li>For the tutorial-relevant fixed points <code>Γ</code>, <code>X</code>, and <code>L</code>, and only for primary subgroup directions (<code>P...</code>) with a microscopic vector distortion on the occupied Wyckoff sites, reconstruct a representative distorted structure directly from <code>DISPLAY DISTORTION</code>.</li>
      <li>Verify the reconstructed structure with a symmetry finder instead of assuming the subgroup table is reproduced automatically.</li>
      <li>Generate diffraction comparisons only for structures that were actually reconstructed.</li>
    </ol>
    """

    limitations = """
    <ul>
      <li>Parameterized k-lines, k-planes, and general k-vectors are documented but not automatically modeled, because a rigorous calculation requires a user-selected k-vector value instead of an arbitrary placeholder.</li>
      <li>The full subgroup catalog is preserved, but the expensive structure/diffraction modeling is limited to the tutorial-relevant <code>Γ</code>, <code>X</code>, and <code>L</code> branches and to the primary <code>P</code>-type directions. Lower-symmetry <code>C</code>, <code>S</code>, generic multidimensional directions, and the <code>W</code>-point family remain documented in the catalog instead of being overinterpreted automatically.</li>
      <li>Some order-parameter directions do not produce an atomic vector distortion for the occupied Wyckoff sites, or they require a more careful supercell embedding than the representative parent-cell reconstruction used here. Those cases are reported explicitly instead of being treated as valid diffraction models.</li>
      <li>Domain counts are reported from ISOTROPY, but the reconstructed diffraction model uses one representative domain rather than all orientational variants.</li>
    </ul>
    """

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Isotropy Scientific Report</title>
  <style>
    :root {{
      --bg: #f4f1ea;
      --panel: #fffdf8;
      --ink: #1f1f1f;
      --muted: #6c6760;
      --accent: #7e2f2f;
      --line: #d7d1c7;
      --ok: #1f6b42;
      --warn: #9a5a00;
      --neutral: #5b6470;
    }}
    body {{
      margin: 0;
      font-family: Georgia, "Times New Roman", serif;
      color: var(--ink);
      background: linear-gradient(180deg, #efe9de 0%, #f8f5ef 22%, #f4f1ea 100%);
      line-height: 1.55;
    }}
    .page {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 32px 24px 64px;
    }}
    h1, h2, h3, h4 {{
      font-family: "Palatino Linotype", "Book Antiqua", Palatino, serif;
      margin: 0 0 12px;
      color: #251d18;
    }}
    h1 {{
      font-size: 2.35rem;
      letter-spacing: 0.01em;
    }}
    h2 {{
      margin-top: 32px;
      padding-top: 16px;
      border-top: 1px solid var(--line);
    }}
    p.lead {{
      max-width: 920px;
      font-size: 1.08rem;
      color: #2b2622;
    }}
    .intro {{
      background: rgba(255,255,255,0.72);
      border: 1px solid rgba(126,47,47,0.18);
      border-radius: 18px;
      padding: 24px;
      box-shadow: 0 12px 28px rgba(58, 43, 32, 0.07);
    }}
    .two-col {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 20px;
      align-items: start;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: var(--panel);
      border: 1px solid var(--line);
      font-size: 0.96rem;
    }}
    th, td {{
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      vertical-align: top;
    }}
    th {{
      background: #f0e8dc;
      text-align: left;
      font-weight: 600;
    }}
    .mode-card {{
      margin: 18px 0;
      background: rgba(255,255,255,0.85);
      border: 1px solid var(--line);
      border-radius: 14px;
      overflow: hidden;
    }}
    .mode-card summary {{
      display: flex;
      justify-content: space-between;
      gap: 16px;
      cursor: pointer;
      list-style: none;
      padding: 18px 20px;
      font-weight: 700;
      background: #f7f2ea;
    }}
    .mode-card summary::-webkit-details-marker {{
      display: none;
    }}
    .mode-card > div,
    .mode-card > h4,
    .mode-card > ul {{
      padding-left: 20px;
      padding-right: 20px;
    }}
    .mode-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 16px;
      padding-top: 18px;
    }}
    .status-ok, .status-warn, .status-neutral {{
      display: inline-block;
      padding: 4px 10px;
      border-radius: 999px;
      font-size: 0.85rem;
      font-weight: 700;
    }}
    .status-ok {{
      background: rgba(31, 107, 66, 0.12);
      color: var(--ok);
    }}
    .status-warn {{
      background: rgba(154, 90, 0, 0.12);
      color: var(--warn);
    }}
    .status-neutral {{
      background: rgba(91, 100, 112, 0.12);
      color: var(--neutral);
    }}
    .figure-block {{
      padding: 12px 20px 0;
    }}
    .figure-block img {{
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: #fff;
    }}
    code {{
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", monospace;
      font-size: 0.92em;
      background: rgba(126, 47, 47, 0.08);
      padding: 0.1em 0.3em;
      border-radius: 4px;
    }}
    @media (max-width: 720px) {{
      .page {{
        padding: 18px 14px 44px;
      }}
      h1 {{
        font-size: 1.8rem;
      }}
    }}
  </style>
</head>
<body>
  <div class="page">
    <section class="intro">
      <h1>ISOTROPY Symmetry and Diffraction Report</h1>
      <p class="lead">
        This report separates what is known directly from the ISOTROPY database from what was actually reconstructed and symmetry-checked from the input CIF.
        That distinction matters physically: an isotropy subgroup is a group-theoretical possibility, while a distorted structure used for diffraction must also be
        embedded consistently on the occupied Wyckoff sites of the parent structure.
      </p>
      {make_html_table(["Metric", "Value"], summary_rows)}
    </section>

    <h2>Parent Structure</h2>
    <div class="two-col">
      <div>{make_html_table(["Site", "Species", "Wyckoff", "Fractional coordinates"], site_rows)}</div>
      <div>
        <h3>Methodology</h3>
        {methodology_notes}
        <h3>Limitations</h3>
        {limitations}
      </div>
    </div>

    <h2>Reciprocal-Space Survey</h2>
    <div class="two-col">
      <div>
        <h3>Fixed k-Points Modeled</h3>
        {make_html_table(["k-point", "Coordinates"], fixed_rows)}
      </div>
      <div>
        <h3>Parameterized k-Manifolds Documented Only</h3>
        {make_html_table(["k-point", "Coordinates"], parameterized_rows)}
      </div>
    </div>

    <h2>Group-Theory Catalog</h2>
    <p>
      Each row below is taken from <code>DISPLAY ISOTROPY</code> for a fixed k-point irrep. The catalog is rigorous as a symmetry listing.
      The diffraction modeling later in this report is limited to cases where <code>DISPLAY DISTORTION</code> yields a microscopic vector mode that can be reconstructed and checked.
    </p>
    {make_html_table(
        ["k-point", "Irrep", "Dir", "Subgroup", "Order parameter", "Basis", "Origin", "Domains"],
        catalog_rows,
    )}

    <h2>Modeled Distortions</h2>
    <p>
      Only the primary <code>P</code>-type subgroup directions on the <code>Γ</code>, <code>X</code>, and <code>L</code> branches were carried through to structure generation, and only when they produced an explicit vector distortion on the occupied Wyckoff sites.
      A status of <code>verified</code> means the reconstructed distorted structure was independently identified with the predicted subgroup symmetry.
    </p>
    {''.join(modeled_sections) or '<p>No microscopic vector distortions were reconstructed.</p>'}
  </div>
</body>
</html>
"""

    with open(REPORT_HTML_PATH, "w", encoding="utf-8") as handle:
        handle.write(html)

    markdown_lines = [
        f"# ISOTROPY Scientific Report: {parent['formula']} ({parent['sg_num']} {parent['sg_symbol']})",
        "",
        "## Summary",
        f"- Fixed k-points analyzed: {report_data['summary']['fixed_kpoint_count']}",
        f"- Parameterized k-manifolds documented only: {report_data['summary']['parameterized_kpoint_count']}",
        f"- Subgroup directions cataloged: {report_data['summary']['catalog_count']}",
        f"- Microscopic vector modes reconstructed: {report_data['summary']['mode_count']}",
        f"- Symmetry-verified distorted structures: {report_data['summary']['verified_mode_count']}",
        "",
        "The detailed version of this report is easier to read in HTML:",
        f"- `scientific_report.html`",
        "",
        "## Key Corrections Relative To The Previous Script",
        "- The workflow now uses `PAGE 1000`, because this ISOTROPY build rejects values larger than 1000.",
        "- Parameterized k-lines, planes, and general points are no longer assigned arbitrary placeholder `KVALUE` values and treated as if they were rigorously analyzed.",
        "- Microscopic distortions are reconstructed directly from `DISPLAY DISTORTION` instead of inventing `x`, `y`, `z` coordinate values.",
        "- The report now distinguishes group-theoretical subgroup listings from symmetry-verified distorted structures.",
        "",
        "## Generated Files",
        "- `scientific_report.html`: readable detailed report",
        "- `scientific_report.json`: structured machine-readable data for the run",
        "- `diffraction_plots/*.png`: diffraction overlays for reconstructed modes",
        "- `subgroups/*.cif`: reconstructed distorted structures",
    ]
    with open(REPORT_PATH, "w", encoding="utf-8") as handle:
        handle.write("\n".join(markdown_lines) + "\n")


def main() -> None:
    parent_cif = resolve_input_path("nacl_standard.cif")
    if len(sys.argv) > 1:
        parent_cif = resolve_input_path(sys.argv[1])

    parent_info, parent_structure = get_parent_info(parent_cif)
    kpoints = get_kpoints(parent_info.sg_num)
    fixed_kpoints = [kp for kp in kpoints if kp.is_fixed]

    catalog = []
    mode_results: list[ModeResult] = []
    seen_signatures: set[str] = set()

    for kpoint in fixed_kpoints:
        irreps = get_irreps(parent_info.sg_num, kpoint.label)
        for irrep in irreps:
            subgroups = get_subgroups(parent_info.sg_num, kpoint.label, irrep)
            for subgroup in subgroups:
                domain_count = get_domain_count(parent_info.sg_num, kpoint.label, irrep, subgroup.direction)
                catalog.append(
                    {
                        "kpoint": kpoint.label,
                        "irrep": irrep,
                        "direction": subgroup.direction,
                        "sg_num": subgroup.sg_num,
                        "symbol": subgroup.symbol,
                        "order_parameter": subgroup.order_parameter,
                        "basis_text": subgroup.basis_text,
                        "origin_text": subgroup.origin_text,
                        "domain_count": domain_count,
                    }
                )

                if kpoint.label not in MODELED_KPOINTS:
                    continue
                if not subgroup.direction.startswith("P"):
                    continue

                notes: list[str] = []
                rows = get_vector_distortion_rows(
                    parent_info.sg_num,
                    kpoint.label,
                    irrep,
                    subgroup.direction,
                    parent_info.wyckoff_letters,
                    domain=None,
                )
                if not rows:
                    notes.append(
                        "No microscopic vector distortion was reported on the occupied Wyckoff sites in the parent conventional cell."
                    )
                    mode_results.append(
                        ModeResult(
                            kpoint=kpoint.label,
                            irrep=irrep,
                            subgroup=subgroup,
                            domain_count=domain_count,
                            displacement_rows=[],
                            displaced_site_count=0,
                            displaced_species=[],
                            max_displacement_norm=0.0,
                            rms_displacement_norm=0.0,
                            verification_status="no_vector_mode",
                            verified_sg_num=None,
                            verified_sg_symbol=None,
                            xrd_plot_file=None,
                            structure_file=None,
                            xrd_summary=None,
                            notes=notes,
                        )
                    )
                    continue

                try:
                    distorted_structure, applied_rows = build_distorted_structure(
                        parent_structure,
                        parent_info,
                        rows,
                        DISTORTION_AMPLITUDE,
                    )
                except ValueError as exc:
                    notes.append(str(exc))
                    mode_results.append(
                        ModeResult(
                            kpoint=kpoint.label,
                            irrep=irrep,
                            subgroup=subgroup,
                            domain_count=domain_count,
                            displacement_rows=[],
                            displaced_site_count=0,
                            displaced_species=[],
                            max_displacement_norm=0.0,
                            rms_displacement_norm=0.0,
                            verification_status="construction_failed",
                            verified_sg_num=None,
                            verified_sg_symbol=None,
                            xrd_plot_file=None,
                            structure_file=None,
                            xrd_summary=None,
                            notes=notes,
                        )
                    )
                    continue

                displaced_site_count, displaced_species, max_norm, rms_norm = summarize_displacements(applied_rows)
                verified_sg_num, verified_sg_symbol = verify_structure(distorted_structure)
                status = "verified" if verified_sg_num == subgroup.sg_num else "mismatch"
                if status == "verified":
                    notes.append("The reconstructed structure reproduces the subgroup symmetry predicted by ISOTROPY.")
                else:
                    notes.append(
                        f"Group-theory prediction is {subgroup.sg_num} {subgroup.symbol}, but the reconstructed parent-cell distortion verifies as {verified_sg_num} {verified_sg_symbol}."
                    )
                    notes.append(
                        "This usually means the subgroup direction requires a more careful supercell embedding, coupled modes, or additional strain/origin handling than the present representative reconstruction."
                    )

                label = f"{kpoint.label}_{irrep}_{subgroup.direction}"
                signature = structure_signature(distorted_structure)
                structure_file = None
                xrd_plot_file = None
                xrd_summary = None
                if status == "verified":
                    if signature in seen_signatures:
                        notes.append(
                            "This verified structure is symmetry-equivalent to a previously saved phase, so no duplicate CIF or diffraction figure was written."
                        )
                    else:
                        seen_signatures.add(signature)
                        structure_file = write_structure_file(distorted_structure, label)
                        xrd_plot_file, xrd_summary = make_xrd_plot(parent_structure, distorted_structure, label)
                else:
                    notes.append(
                        "Because the representative reconstruction did not verify the predicted subgroup symmetry, no diffraction figure is claimed for this mode."
                    )
                mode_results.append(
                    ModeResult(
                        kpoint=kpoint.label,
                        irrep=irrep,
                        subgroup=subgroup,
                        domain_count=domain_count,
                        displacement_rows=applied_rows,
                        displaced_site_count=displaced_site_count,
                        displaced_species=displaced_species,
                        max_displacement_norm=max_norm,
                        rms_displacement_norm=rms_norm,
                        verification_status=status,
                        verified_sg_num=verified_sg_num,
                        verified_sg_symbol=verified_sg_symbol,
                        xrd_plot_file=xrd_plot_file,
                        structure_file=structure_file,
                        xrd_summary=xrd_summary,
                        notes=notes,
                    )
                )

    report_data = {
        "parent": {
            **asdict(parent_info),
            "site_summaries": [asdict(site) for site in parent_info.site_summaries],
        },
        "kpoints": [asdict(kpoint) for kpoint in kpoints],
        "catalog": catalog,
        "mode_results": [
            {
                **asdict(result),
                "subgroup": asdict(result.subgroup),
                "displacement_rows": [asdict(row) for row in result.displacement_rows],
            }
            for result in mode_results
        ],
        "summary": {
            "fixed_kpoint_count": len(fixed_kpoints),
            "parameterized_kpoint_count": len([kp for kp in kpoints if not kp.is_fixed]),
            "catalog_count": len(catalog),
            "mode_count": len(mode_results),
            "verified_mode_count": len([result for result in mode_results if result.verification_status == "verified"]),
        },
    }
    generate_reports(report_data)
    print(f"HTML report: {REPORT_HTML_PATH}")
    print(f"Markdown summary: {REPORT_PATH}")
    print(f"Structured data: {REPORT_JSON_PATH}")


if __name__ == "__main__":
    main()

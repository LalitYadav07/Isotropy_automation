from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import asdict
from pathlib import Path

try:
    from pymatgen.io.cif import CifParser
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
except ImportError:  # pragma: no cover - optional dependency for verification only
    CifParser = None
    SpacegroupAnalyzer = None

try:
    from _isotropy_paths import TUTORIAL_REPORTS_DIR, TUTORIALS_DIR, configure_environment, resolve_input_path
    from isodistort_distortion_parser import DistortionReportData, parse_distortion_file
except ModuleNotFoundError:
    from isotropy_project.scripts._isotropy_paths import (
        TUTORIAL_REPORTS_DIR,
        TUTORIALS_DIR,
        configure_environment,
        resolve_input_path,
    )
    from isotropy_project.scripts.isodistort_distortion_parser import DistortionReportData, parse_distortion_file


configure_environment()


TUTORIAL_METADATA = {
    "lamno3-magchild": {
        "exercise": "Exercise 3 + Exercise 8",
        "headline": "Mode decomposition and quantitative amplitudes in antiferromagnetic LaMnO3",
        "summary": (
            "This example is not a search from scratch. It is a symmetry-mode decomposition of a known "
            "commensurate magnetic child structure relative to a cubic parent."
        ),
        "objectives": [
            "Understand how a known child structure is expanded in the symmetry-mode basis of the parent.",
            "Separate primary magnetic order from secondary octahedral rotations, Jahn-Teller distortions, and strain.",
            "Interpret ISODISTORT amplitudes as root-summed-squared local order parameters over the full supercell.",
        ],
        "workflow": [
            "Import `lamno3-cubic.cif` as the parent and enable Mn magnetic moments.",
            "Use Method 4 with `lamno3-magchild.cif` as the child structure.",
            "Keep the guessed basis `{(1,0,-1),(0,2,0),(1,0,1)}` and automatic origin detection.",
            "Choose nearest-neighbor atom matching, inspect the mode decomposition, and save the distortion file.",
        ],
        "interpretation": (
            "The tutorial emphasizes that the experimental magnetic structure is carried by the `mX5+` order parameter, "
            "while multiple displacive irreps and strain modes remain symmetry-allowed and acquire nonzero amplitudes in the decomposition. "
            "This is the right example for learning how ISODISTORT defines mode amplitudes quantitatively."
        ),
    },
    "dymn6ge6": {
        "exercise": "Exercise 5",
        "headline": "Incommensurate magnetic DyMn6Ge6 and superspace verification",
        "summary": (
            "This example superposes a commensurate ferromagnetic irrep with an incommensurate magnetic irrep to form a conical magnetic structure."
        ),
        "objectives": [
            "Build an incommensurate magnetic structure with a superspace-group description.",
            "Track the difference between MSSG, BMSG, BFSG, and FSSG viewpoints.",
            "Use exported CIF data as a verifiable bridge to FINDSSG and ISOCIF workflows.",
        ],
        "workflow": [
            "Import `dymn6ge6-parent.cif` and disable displacements while enabling Dy and Mn magnetic moments.",
            "Use Method 2 with two irreps: `mGM2+` at `GM` and `mDT6` at `DT(0,0,g)` with `g = 0.165`.",
            "Choose OPD `(a|b,0,0,0)` to obtain MSSG `177.1.24.2.m153.1 P62'2'(0,0,g)h00`.",
            "Enter the tutorial amplitudes, inspect the mode details, and export the modulated CIF.",
        ],
        "interpretation": (
            "The file encodes both the zero-wave-vector ferromagnetic contribution and the incommensurate conical magnetic modulation. "
            "It is a strong test case for any report because conventional space-group verification is insufficient; the superspace description is essential."
        ),
    },
    "tbmno3": {
        "exercise": "Exercise 6",
        "headline": "Cycloidal magnetic order in TbMnO3 before explicit polar displacements",
        "summary": (
            "This distortion file captures the primary incommensurate magnetic OPD that produces the multiferroic isotropy subgroup."
        ),
        "objectives": [
            "Understand how two incommensurate magnetic irreps combine into a cycloid.",
            "See that secondary ferroelectric order is symmetry-allowed even before explicit displacive amplitudes are entered.",
            "Connect magnetic OPD selection to the loss of inversion symmetry and multiferroic behavior.",
        ],
        "workflow": [
            "Import `tbmno3-pbnm.cif`, keep strain and displacements enabled, and enable Mn magnetic moments.",
            "Use Method 2 with two `SM(0,a,0)` irreps at `a = 0.27`: `mSM3` and `mSM2`.",
            "Choose OPD `(a,0|b,0)` in the `ba-c` orthorhombic setting with basic superspace setting.",
            "Inspect the magnetic mode content and export the distortion before or after adding polar displacements.",
        ],
        "interpretation": (
            "This file is useful for separating symmetry from numerics. The chosen OPD already fixes a polar superspace subgroup, "
            "but the stored coefficients show only the primary magnetic modes. The ferroelectric displacive response is symmetry-allowed but not yet populated numerically."
        ),
    },
    "tbmno3-distortion2": {
        "exercise": "Exercise 6 extension",
        "headline": "Cycloidal TbMnO3 with explicit secondary ferroelectric displacements",
        "summary": (
            "This second TbMnO3 file appears to be the same magnetic solution after adding explicit displacive amplitudes for the secondary polar response."
        ),
        "objectives": [
            "Contrast a purely magnetic OPD selection with a distortion file that also stores secondary displacive amplitudes.",
            "Show how ISODISTORT keeps primary and secondary order parameters inside one symmetry-consistent description.",
            "Document the difference between subgroup symmetry and a user-populated numerical realization of that subgroup.",
        ],
        "workflow": [
            "Repeat the main TbMnO3 exercise through the cycloidal `mSM3 + mSM2` OPD.",
            "In ISOVIZ, add small `GM4-` polar displacements of cations against oxygen.",
            "Transfer those amplitudes back to the distortion page and save a second distortion file.",
        ],
        "interpretation": (
            "The presence of nonzero displacive coefficients makes this file closer to a physically populated multiferroic state than the first TbMnO3 file. "
            "The comparison between the two files is itself pedagogically useful."
        ),
    },
    "skyrmion": {
        "exercise": "Exercise 7",
        "headline": "Multi-k superspace skyrmion lattice on a Fe monolayer",
        "summary": (
            "This is the most symmetry-rich example in the archive: a two-modulation, three-arm, multi-k magnetic superspace construction."
        ),
        "objectives": [
            "Work with an OPD that keeps symmetry-equivalent modulation waves locked together.",
            "Understand why a three-arm star of k can still be described in a `(3+2)`D superspace.",
            "See how occupational, displacive, magnetic, and strain channels can all coexist in one symmetry scaffold.",
        ],
        "workflow": [
            "Import `skyrmion-parent.cif`, enable Fe occupational and magnetic order parameters.",
            "Use Method 2 with `LD(a,a,0)`, `a = 0.1`, and `d = 2` independent modulations.",
            "Choose `mLD3` or the `mLD3 + mLD4` superposition depending on whether you want a vortex lattice or a full skyrmion texture.",
            "Save the distortion file after setting the magnetic amplitudes of interest.",
        ],
        "interpretation": (
            "The saved file in this archive is especially informative when its coefficients are zero, because it still captures the full symmetry scaffold: "
            "multi-arm k-star, superspace subgroup, and all allowed secondary channels. That is useful for documentation even before a numerical state is populated."
        ),
    },
}


def html_escape(text: object) -> str:
    value = str(text)
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def relpath(from_dir: Path, to_path: Path) -> str:
    return os.path.relpath(to_path, start=from_dir)


def parse_order_parameter_details(order_parameter: str) -> dict[str, str]:
    details: dict[str, str] = {"raw": order_parameter}
    lead = order_parameter.split(", basis=", 1)[0].strip()
    basis_match = re.search(r"basis=\{([^}]*)\}", order_parameter)
    origin_match = re.search(r"origin=\(([^)]*)\)", order_parameter)
    size_match = re.search(r"\bs=([^,]+)", order_parameter)
    index_match = re.search(r"\bi=([^,]+)", order_parameter)
    k_match = re.search(r"k-active=\s*(.+)$", order_parameter)

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
    if k_match:
        details["k_active"] = k_match.group(1).strip()
    return details


def distortion_type_summary(data: DistortionReportData) -> list[str]:
    labels: list[str] = []
    if "displacive" in data.classification:
        labels.append("atomic displacements")
    if "magnetic" in data.classification:
        labels.append("magnetic moments")
    if "occupational" in data.classification:
        labels.append("occupational order")
    if "rotational" in data.classification:
        labels.append("rotational modes")
    if "strain" in data.classification:
        labels.append("strain")
    return labels


def active_mode_channels(data: DistortionReportData) -> list[dict[str, object]]:
    channels = []
    for family, summary in data.mode_summaries.items():
        if summary["nonzero_mode_count"] == 0:
            continue
        channels.append(
            {
                "family": family,
                "nonzero_mode_count": summary["nonzero_mode_count"],
                "rss": summary["rss"],
                "max_abs": summary["max_abs"],
            }
        )
    channels.sort(key=lambda item: item["rss"], reverse=True)
    return channels


def verify_cif_artifacts(source_path: Path, associated_files: list[str]) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for filename in associated_files:
        path = source_path.parent / filename
        entry: dict[str, object] = {
            "filename": filename,
            "exists": path.exists(),
            "suffix": path.suffix.lower(),
        }
        if not path.exists() or path.suffix.lower() != ".cif":
            results.append(entry)
            continue
        if CifParser is None or SpacegroupAnalyzer is None:
            entry["parse_status"] = "verification dependency unavailable"
            results.append(entry)
            continue
        try:
            structure = CifParser(str(path)).parse_structures(primitive=False)[0]
            analyzer = SpacegroupAnalyzer(structure, symprec=1e-2)
            entry["parse_status"] = "parsed"
            entry["formula"] = structure.composition.reduced_formula
            entry["space_group"] = f"{analyzer.get_space_group_number()} {analyzer.get_space_group_symbol()}"
        except Exception as exc:  # pragma: no cover - depends on parser support
            entry["parse_status"] = f"not parsed: {exc.__class__.__name__}"
        results.append(entry)
    return results


def verification_summary(data: DistortionReportData, source_path: Path) -> dict[str, object]:
    cif_results = verify_cif_artifacts(source_path, data.associated_files)
    parsed_cifs = [item for item in cif_results if item.get("parse_status") == "parsed"]
    failed_cifs = [
        item
        for item in cif_results
        if item.get("suffix") == ".cif" and item.get("parse_status") not in {None, "parsed"}
    ]
    return {
        "internal_checks": data.internal_checks,
        "cif_results": cif_results,
        "parsed_cif_count": len(parsed_cifs),
        "failed_cif_count": len(failed_cifs),
    }


def build_example_record(data: DistortionReportData, source_path: Path, output_dir: Path) -> dict[str, object]:
    metadata = TUTORIAL_METADATA.get(data.slug, {})
    order_details = parse_order_parameter_details(data.order_parameter)
    verification = verification_summary(data, source_path)
    source_rel = relpath(output_dir, source_path)
    associated_links = [
        {"filename": filename, "href": relpath(output_dir, source_path.parent / filename)}
        for filename in data.associated_files
    ]

    channel_rows = []
    for family, summary in data.mode_summaries.items():
        channel_rows.append(
            {
                "family": family,
                "total_mode_count": summary["total_mode_count"],
                "nonzero_mode_count": summary["nonzero_mode_count"],
                "rss": summary["rss"],
                "max_abs": summary["max_abs"],
                "consistency_ok": summary["consistency_ok"],
            }
        )

    return {
        "slug": data.slug,
        "title": metadata.get("headline", data.title.title()),
        "exercise": metadata.get("exercise", "Tutorial example"),
        "summary": metadata.get("summary", ""),
        "objectives": metadata.get("objectives", []),
        "workflow": metadata.get("workflow", []),
        "interpretation": metadata.get("interpretation", ""),
        "source_file": source_rel,
        "source_filename": source_path.name,
        "classification": data.classification,
        "parent": {
            "label": data.parent_string,
            "formula_estimate": data.parent_formula_estimate,
            "setting": data.parent_setting,
            "lattice": data.lattice_string,
            "wyckoff_sites": data.wyckoff_strings,
        },
        "distortion_types": distortion_type_summary(data),
        "primary_selection": {
            "k_vectors": data.k_vectors,
            "irreps": data.irreps,
            "details": order_details,
        },
        "mode_overview": {
            "mode_counts": data.mode_counts,
            "channels": channel_rows,
            "active_channels": active_mode_channels(data),
            "primary_irreps": data.primary_irreps,
        },
        "child_basic_structure": data.child_basic_structure,
        "artifacts": associated_links,
        "verification": verification,
        "raw_data": asdict(data),
    }


def render_badges(items: list[str]) -> str:
    return "".join(f'<span class="badge">{html_escape(item)}</span>' for item in items)


def render_table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html_escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def render_list(items: list[str]) -> str:
    if not items:
        return "<p class=\"muted\">None recorded.</p>"
    return "<ul>" + "".join(f"<li>{html_escape(item)}</li>" for item in items) + "</ul>"


def render_cif_results(results: list[dict[str, object]]) -> str:
    rows = []
    for item in results:
        rows.append(
            [
                html_escape(item["filename"]),
                html_escape(item.get("parse_status", "not checked")),
                html_escape(item.get("formula", "")),
                html_escape(item.get("space_group", "")),
            ]
        )
    return render_table(["Artifact", "Parse status", "Formula", "Average-structure SG"], rows or [["None", "", "", ""]])


def render_primary_irreps(primary_irreps: list[dict[str, object]]) -> str:
    if not primary_irreps:
        return "<p class=\"muted\">No explicit primary-irrep list is stored in this distortion file. That is expected for mode-decomposition exports.</p>"
    rows = []
    for item in primary_irreps:
        categories = []
        for family, summary in item.get("category_summaries", {}).items():
            if summary["nonzero_count"] > 0:
                categories.append(f"{family}: rss={summary['rss']:.3f}")
        rows.append(
            [
                str(item["index"]),
                html_escape(item.get("irrep", "")),
                html_escape(item.get("k_vector", "")),
                html_escape(item.get("irrep_number", "")),
                html_escape(", ".join(categories) or "no nonzero stored coefficients"),
            ]
        )
    return render_table(["#", "Irrep", "k vector", "Internal irrep number", "Stored amplitude summary"], rows)


def render_channel_table(channels: list[dict[str, object]]) -> str:
    rows = []
    for item in channels:
        rows.append(
            [
                html_escape(item["family"]),
                str(item["total_mode_count"]),
                str(item["nonzero_mode_count"]),
                f"{item['rss']:.4f}",
                f"{item['max_abs']:.4f}",
                "yes" if item["consistency_ok"] else "no",
            ]
        )
    return render_table(
        ["Channel", "Total modes", "Active coefficients", "RSS coefficient", "Max |coef|", "Length checks"],
        rows,
    )


def render_associated_artifacts(artifacts: list[dict[str, str]]) -> str:
    if not artifacts:
        return "<p class=\"muted\">No sibling tutorial artifacts were detected.</p>"
    items = []
    for artifact in artifacts:
        items.append(
            f'<li><a href="{html_escape(artifact["href"])}">{html_escape(artifact["filename"])}</a></li>'
        )
    return "<ul>" + "".join(items) + "</ul>"


def page_css() -> str:
    return """
    :root {
      --bg: #f6f3ec;
      --panel: #fffdf8;
      --ink: #1f1d1a;
      --muted: #6d655c;
      --accent: #0a5c57;
      --accent-soft: #d8ebe6;
      --line: #ddd3c7;
      --warn: #8a5a14;
      --good: #256b34;
      --mono: "IBM Plex Mono", "SFMono-Regular", Consolas, monospace;
      --serif: "Iowan Old Style", "Palatino Linotype", "Book Antiqua", Palatino, serif;
      --sans: "Avenir Next", "Segoe UI", sans-serif;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      background:
        radial-gradient(circle at top left, rgba(10,92,87,0.10), transparent 22rem),
        linear-gradient(180deg, #faf8f2 0%, var(--bg) 100%);
      color: var(--ink);
      font-family: var(--sans);
      line-height: 1.55;
    }
    main { max-width: 1180px; margin: 0 auto; padding: 2rem 1.25rem 4rem; }
    h1, h2, h3, h4 { font-family: var(--serif); line-height: 1.15; margin: 0 0 0.75rem; }
    h1 { font-size: clamp(2.2rem, 4vw, 3.5rem); }
    h2 { font-size: 1.7rem; margin-top: 2rem; }
    h3 { font-size: 1.25rem; margin-top: 1.5rem; }
    p, li { font-size: 1rem; }
    .lead { font-size: 1.1rem; max-width: 75ch; color: var(--muted); }
    .panel {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 1rem 1.1rem;
      box-shadow: 0 12px 30px rgba(48, 36, 25, 0.05);
    }
    .grid { display: grid; gap: 1rem; }
    .grid.two { grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); }
    .grid.three { grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); }
    .card-link { color: inherit; text-decoration: none; }
    .badge {
      display: inline-block;
      margin: 0 0.45rem 0.45rem 0;
      padding: 0.22rem 0.6rem;
      border-radius: 999px;
      background: var(--accent-soft);
      color: var(--accent);
      font-size: 0.82rem;
      font-weight: 600;
    }
    .metric {
      font-family: var(--mono);
      font-size: 0.92rem;
      color: var(--muted);
    }
    .muted { color: var(--muted); }
    code, .formula { font-family: var(--mono); font-size: 0.92em; }
    table {
      width: 100%;
      border-collapse: collapse;
      font-size: 0.95rem;
      background: white;
      border: 1px solid var(--line);
    }
    th, td {
      padding: 0.65rem 0.6rem;
      border-bottom: 1px solid var(--line);
      vertical-align: top;
      text-align: left;
    }
    th { background: #f4efe6; font-weight: 700; }
    .hero {
      padding: 1.3rem 1.4rem;
      border-radius: 24px;
      background: linear-gradient(135deg, rgba(10,92,87,0.12), rgba(138,90,20,0.08));
      border: 1px solid rgba(10,92,87,0.18);
      margin-bottom: 1.5rem;
    }
    .section-note {
      border-left: 4px solid var(--accent);
      padding-left: 0.9rem;
      color: var(--muted);
      max-width: 75ch;
    }
    .status-good { color: var(--good); font-weight: 700; }
    .status-warn { color: var(--warn); font-weight: 700; }
    .topbar {
      display: flex;
      justify-content: space-between;
      gap: 1rem;
      align-items: flex-start;
      flex-wrap: wrap;
      margin-bottom: 1rem;
    }
    a { color: var(--accent); }
    @media (max-width: 720px) {
      main { padding: 1rem 0.9rem 3rem; }
      .panel { padding: 0.9rem; }
    }
    """


def render_example_html(record: dict[str, object], output_path: Path) -> str:
    verification = record["verification"]
    checks = verification["internal_checks"]
    order_details = record["primary_selection"]["details"]
    parent = record["parent"]
    child = record["child_basic_structure"] or {}
    check_rows = [
        ["Associated files present", "yes" if checks["associated_files_all_present"] else "no"],
        ["Coefficient tables internally consistent", "yes" if checks["coefficient_lengths_consistent"] else "no"],
        ["Primary-irrep alignment inferred", "yes" if checks["primary_irrep_alignment"] else "no"],
        ["Atoms section embedded", "yes" if checks["atoms_section_present"] else "no"],
    ]
    if checks["missing_associated_files"]:
        check_rows.append(["Missing associated files", html_escape(", ".join(checks["missing_associated_files"]))])

    child_rows = []
    if child:
        child_rows = [
            ["Basic-cell atom records", str(child.get("atom_count", 0))],
            [
                "Basic-cell lattice",
                ", ".join(f"{value:.5f}" for value in child.get("lattice_parameters", [])),
            ],
        ]

    source_link = html_escape(record["source_file"])
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html_escape(record["title"])}</title>
  <style>{page_css()}</style>
</head>
<body>
  <main>
    <div class="topbar">
      <div>
        <p class="metric"><a href="index.html">Tutorial index</a></p>
        <h1>{html_escape(record["title"])}</h1>
        <p class="lead">{html_escape(record["summary"])}</p>
      </div>
      <div class="panel">
        <div class="metric">{html_escape(record["exercise"])}</div>
        <div>{render_badges(record["classification"])}</div>
        <p class="metric">Source: <a href="{source_link}">{html_escape(record["source_filename"])}</a></p>
      </div>
    </div>

    <section class="hero">
      <h2>What This Example Teaches</h2>
      {render_list(record["objectives"])}
    </section>

    <section class="grid two">
      <div class="panel">
        <h2>Parent Reference</h2>
        <p><strong>Parent symmetry:</strong> <span class="formula">{html_escape(parent["label"])}</span></p>
        <p><strong>Estimated parent formula:</strong> <span class="formula">{html_escape(parent["formula_estimate"] or "not reconstructed")}</span></p>
        <p><strong>Parent setting:</strong> {html_escape(parent["setting"] or "default / unspecified")}</p>
        <p><strong>Parent lattice:</strong> <span class="formula">{html_escape(parent["lattice"])}</span></p>
        <h3>Wyckoff Sites Used In The Tutorial</h3>
        {render_list(parent["wyckoff_sites"])}
      </div>

      <div class="panel">
        <h2>OPD And Subgroup</h2>
        <p><strong>Stored selection:</strong> <span class="formula">{html_escape(order_details.get("subgroup_selection", order_details["raw"]))}</span></p>
        <p><strong>OPD label:</strong> <span class="formula">{html_escape(order_details.get("opd_label", "not explicit"))}</span></p>
        <p><strong>OPD vector:</strong> <span class="formula">{html_escape(order_details.get("opd_vector", "not explicit"))}</span></p>
        <p><strong>Basis:</strong> <span class="formula">{html_escape(order_details.get("basis", "not embedded"))}</span></p>
        <p><strong>Origin:</strong> <span class="formula">{html_escape(order_details.get("origin", "not embedded"))}</span></p>
        <p><strong>Size index s:</strong> <span class="formula">{html_escape(order_details.get("size_index", "not embedded"))}</span></p>
        <p><strong>Symmetry index i:</strong> <span class="formula">{html_escape(order_details.get("symmetry_index", "not embedded"))}</span></p>
        <p><strong>Active k vectors in OPD:</strong> <span class="formula">{html_escape(order_details.get("k_active", "not embedded"))}</span></p>
      </div>
    </section>

    <section class="grid two">
      <div class="panel">
        <h2>Primary Selection</h2>
        <p><strong>Primary k vectors:</strong></p>
        {render_list(record["primary_selection"]["k_vectors"])}
        <p><strong>Primary irreps:</strong></p>
        {render_list(record["primary_selection"]["irreps"])}
        <p class="section-note">{html_escape(record["interpretation"])}</p>
      </div>

      <div class="panel">
        <h2>Stored Distortion Types</h2>
        {render_list(record["distortion_types"])}
        <h3>Tutorial Reproduction Path</h3>
        {render_list(record["workflow"])}
      </div>
    </section>

    <section class="panel">
      <h2>Mode-Amplitude Overview</h2>
      <p class="lead">The table below summarizes the coefficients stored in the distortion file, not a re-derived physical refinement. This distinction matters for zero-filled symmetry scaffolds and for files where only primary modes were populated manually.</p>
      {render_channel_table(record["mode_overview"]["channels"])}
      <h3>Primary-Irrep Breakdown</h3>
      {render_primary_irreps(record["mode_overview"]["primary_irreps"])}
    </section>

    <section class="grid two">
      <div class="panel">
        <h2>Verification Scope</h2>
        <p class="section-note">These checks are intentionally conservative. They verify internal consistency of the exported distortion file and test whether companion CIF artifacts are parseable as average crystal structures. They do not replace magnetic or superspace symmetry validation in FINDSSG, ISOCIF, or JANA-style tools.</p>
        {render_table(["Check", "Result"], check_rows)}
        <h3>Artifact Parsing</h3>
        {render_cif_results(verification["cif_results"])}
      </div>

      <div class="panel">
        <h2>Exported Artifacts</h2>
        {render_associated_artifacts(record["artifacts"])}
        <h3>Embedded Child Structure Snapshot</h3>
        {render_table(["Field", "Value"], child_rows or [["Embedded structure", "not included in this file"]])}
      </div>
    </section>
  </main>
</body>
</html>
"""


def render_index_html(records: list[dict[str, object]]) -> str:
    cards = []
    for record in records:
        checks = record["verification"]["internal_checks"]
        status_class = "status-good" if checks["coefficient_lengths_consistent"] else "status-warn"
        cards.append(
            f"""
            <a class="card-link" href="{html_escape(record['slug'])}.html">
              <article class="panel">
                <div class="topbar">
                  <div>
                    <div class="metric">{html_escape(record['exercise'])}</div>
                    <h3>{html_escape(record['title'])}</h3>
                  </div>
                  <div class="{status_class}">{'checked' if checks['coefficient_lengths_consistent'] else 'needs review'}</div>
                </div>
                <p>{html_escape(record['summary'])}</p>
                <div>{render_badges(record['classification'])}</div>
              </article>
            </a>
            """
        )

    primer_rows = [
        ["Parent structure", "The higher-symmetry reference structure against which irreps and symmetry modes are defined."],
        ["Irrep", "An irreducible representation of the parent symmetry at a chosen k point or star of k."],
        ["OPD", "Order-parameter direction. It selects how the components of a multi-dimensional irrep are populated."],
        ["Isotropy subgroup", "The child symmetry retained by the selected order parameter, identified by subgroup type plus basis and origin."],
        ["Mode amplitude", "ISODISTORT defines amplitudes as root-summed-squared local order parameters over all affected supercell atoms."],
        ["Superspace group", "The correct symmetry language for incommensurate structures; ordinary 3D space groups are not sufficient."],
    ]

    scope_rows = [
        ["What is verified automatically", "Distortion-file consistency, artifact presence, and CIF parseability when a conventional parser can read the export."],
        ["What is not claimed automatically", "Complete magnetic or superspace symmetry proof, domain enumeration, or a full refinement-quality physical interpretation from CIF parsing alone."],
        ["Why this matters", "Several tutorial examples are magnetic or incommensurate, so a simple space-group finder can under-describe the real symmetry content."],
    ]

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ISODISTORT Tutorial Reports</title>
  <style>{page_css()}</style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>ISODISTORT Tutorial Reports</h1>
      <p class="lead">This report set turns the completed IWMC2024 tutorial artifacts into a self-guided reference. Each page separates three things that were previously conflated: the symmetry idea the tutorial is teaching, the exact data stored in the exported distortion file, and the narrower set of facts that can be independently checked from companion artifacts.</p>
    </section>

    <section class="grid two">
      <div class="panel">
        <h2>How To Read These Reports</h2>
        {render_table(["Concept", "Why it matters"], primer_rows)}
      </div>
      <div class="panel">
        <h2>Verification Policy</h2>
        {render_table(["Scope", "Meaning"], scope_rows)}
      </div>
    </section>

    <section>
      <h2>Examples</h2>
      <div class="grid two">
        {''.join(cards)}
      </div>
    </section>
  </main>
</body>
</html>
"""


def render_example_markdown(record: dict[str, object]) -> str:
    lines = [
        f"# {record['title']}",
        "",
        f"Exercise: {record['exercise']}",
        "",
        record["summary"],
        "",
        "## Objectives",
    ]
    lines.extend(f"- {item}" for item in record["objectives"])
    lines.extend(
        [
            "",
            "## Parent",
            f"- Parent symmetry: `{record['parent']['label']}`",
            f"- Estimated formula: `{record['parent']['formula_estimate'] or 'not reconstructed'}`",
            f"- Parent lattice: `{record['parent']['lattice']}`",
            "",
            "## OPD And Subgroup",
            f"- Stored selection: `{record['primary_selection']['details'].get('subgroup_selection', record['primary_selection']['details']['raw'])}`",
            f"- Basis: `{record['primary_selection']['details'].get('basis', 'not embedded')}`",
            f"- Origin: `{record['primary_selection']['details'].get('origin', 'not embedded')}`",
            f"- Active k vectors: `{record['primary_selection']['details'].get('k_active', 'not embedded')}`",
            "",
            "## Interpretation",
            record["interpretation"],
            "",
            "## Verification",
            f"- Associated files present: `{record['verification']['internal_checks']['associated_files_all_present']}`",
            f"- Coefficient lengths consistent: `{record['verification']['internal_checks']['coefficient_lengths_consistent']}`",
            f"- Parsed CIF count: `{record['verification']['parsed_cif_count']}`",
        ]
    )
    return "\n".join(lines) + "\n"


def load_distortion_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return [input_path]
    return sorted(input_path.glob("*-distortion*.txt"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate tutorial-style ISODISTORT reports.")
    parser.add_argument(
        "--input",
        default=str(TUTORIALS_DIR / "exercises-completed"),
        help="Distortion file or directory containing distortion files.",
    )
    parser.add_argument(
        "--output",
        default=str(TUTORIAL_REPORTS_DIR),
        help="Directory for generated tutorial reports.",
    )
    args = parser.parse_args()

    input_path = resolve_input_path(args.input)
    output_dir = Path(args.output)
    if not output_dir.is_absolute():
        output_dir = (Path.cwd() / output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for source_path in load_distortion_files(input_path):
        data = parse_distortion_file(source_path)
        records.append(build_example_record(data, source_path, output_dir))

    records.sort(key=lambda item: (item["exercise"], item["title"]))

    index_json_path = output_dir / "index.json"
    index_md_path = output_dir / "index.md"
    index_html_path = output_dir / "index.html"

    with index_json_path.open("w", encoding="utf-8") as handle:
        json.dump(records, handle, indent=2)

    index_md_lines = [
        "# ISODISTORT Tutorial Reports",
        "",
        "Generated examples:",
    ]
    index_md_lines.extend(f"- [{record['title']}]({record['slug']}.md)" for record in records)
    index_md_path.write_text("\n".join(index_md_lines) + "\n", encoding="utf-8")
    index_html_path.write_text(render_index_html(records), encoding="utf-8")

    for record in records:
        (output_dir / f"{record['slug']}.json").write_text(
            json.dumps(record, indent=2),
            encoding="utf-8",
        )
        (output_dir / f"{record['slug']}.md").write_text(
            render_example_markdown(record),
            encoding="utf-8",
        )
        (output_dir / f"{record['slug']}.html").write_text(
            render_example_html(record, output_dir / f"{record['slug']}.html"),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()

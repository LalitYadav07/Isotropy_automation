# Repository Guidance

## Scope

This repo contains local automation around the BYU ISOTROPY suite for structural phase-transition analysis.

The original package is vendored under `vendor/isobyu/`. Do not modify those files unless the user explicitly asks for vendor changes.

## Main Entry Points

- `scripts/discover_distortion_signatures.py`
  Current main workflow for parent-side discovery from a CIF.
- `scripts/generate_tutorial_reports.py`
  Builds tutorial-style reports from completed ISODISTORT exports.
- `docs/isotropy_transition_model.md`
  Current conceptual model and terminology.
- `README.md`
  Project layout and usage notes.

## What The Current Discovery Workflow Does

The current `discover_distortion_signatures.py` workflow:

1. Standardizes the parent CIF.
2. Enumerates fixed-k single-irrep isotropy subgroups with local `iso`.
3. Screens fixed-k irreps for vector-active structural channels on the occupied parent sites.
4. Builds canonical single-irrep probe structures for primary `P...` directions.
5. Verifies the resulting subgroup geometrically and symmetrically.
6. Simulates coarse powder-diffraction signatures.
7. Optionally catalogs coupled-primary subgroup intersections with `DISPLAY ISOTROPY COUPLED`.
8. Reports transition-class coverage explicitly.

## Important Physics / Symmetry Notes

- `vector_active` means only that symmetry allows a microscopic displacement channel on the occupied parent sites.
- `vector_active` does not imply a soft mode, energetic instability, or experimental relevance.
- `DISPLAY DISTORTION` rows are not always one-row-per-atom. They can be projected-vector contributions in the chosen display cell.
- Those contributions must be transformed into the subgroup basis and merged on the same child-cell site before constructing a probe structure.
- Reconstructive transitions are not discoverable from a parent CIF alone. They require parent-child comparison and likely `comsubs`.
- Incommensurate, one-arm/kernel, and scalar/occupational workflows are not yet automated here.

## Current Practical Gaps

- Parent standardization currently uses `pymatgen` rather than a fully robust direct `findsym` automation path.
- Coupled-primary search is catalog-only at present; coupled child-structure generation is not implemented.
- Parameterized commensurate k-lines / k-planes are documented but not expanded automatically.
- Incommensurate superspace handling is not automated.
- Occupational / scalar structural order is not automated.
- Reconstructive analysis is not yet wrapped in a dedicated script.

## Suggested Next Steps

If continuing the project, the best next development targets are:

1. add a parent+child reconstructive workflow around `comsubs`
2. automate one-arm and parameterized-k searches
3. add scalar / occupational structural channels
4. build coupled-primary child-structure generation and verification

## Setup

Use Python with the dependencies in `requirements.txt`.

Example:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The scripts set `ISODATA` and `PATH` for the vendored executables through `scripts/_isotropy_paths.py`.

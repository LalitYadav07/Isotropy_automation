# Isotropy Project Layout

This folder separates the original BYU ISOTROPY package from the local scripts and generated outputs used to study it.

## Structure

- `vendor/isobyu/`: original ISOTROPY executables, manuals, data tables, and sample files
- `scripts/`: local Python wrappers and debugging utilities
- `inputs/`: example CIFs and hand-written test inputs
- `isotutorials/`: extracted IWMC2024 tutorial materials and completed examples
- `notebooks/`: exploratory notebook work
- `docs/`: local notes on ISOTROPY concepts and workflow assumptions
- `outputs/subgroups/`: generated subgroup CIFs
- `outputs/diffraction_plots/`: generated diffraction plots
- `outputs/logs/`: saved run logs
- `outputs/scientific_report.md`: generated report
- `outputs/tutorial_reports/`: tutorial-style HTML/JSON/Markdown reports generated from ISODISTORT distortion files
- `outputs/discovery_runs/`: per-material subgroup catalogs, probe structures, signature plots, and discovery reports
- `scratch/`: temporary files created while running the scripts

## Script behavior

The scripts in `scripts/` now resolve paths relative to this project root and leave `vendor/isobyu/` unchanged.

Important entry points:

- `scripts/discover_distortion_signatures.py`: parent-side structural distortion discovery
- `scripts/generate_tutorial_reports.py`: distortion-file explanation/report generation
- `scripts/explain_phase_transition.py`: transition-family classifier and workflow router
- `scripts/reconstructive_transition_search.py`: parent+child reconstructive/common-subgroup search with `comsubs`
- `scripts/magnetic_workflows.py`: dedicated magnetic and superspace workflow family entry points

## Tutorial reports

Use `./.venv/bin/python isotropy_project/scripts/generate_tutorial_reports.py` to turn the completed ISODISTORT distortion files in `isotutorials/exercises-completed/` into readable reports.

The generated tutorial reports are designed to:

- explain the parent symmetry, selected irreps, OPD, subgroup basis, and origin
- summarize which distortion channels are actually populated in the saved file
- distinguish between tutorial intent, exported facts, and conservative automated verification
- keep magnetic and incommensurate examples honest by not overstating what a conventional CIF parser can verify

Magnetic and superspace tutorial examples now also have a dedicated regression-fixture set in `fixtures/magnetic_tutorial_regression.json`, documented in `docs/magnetic_validation.md` and enforced by `python -m unittest tests.test_magnetic_tutorial_regression`.

The older `simulate_diffraction.py` workflow remains an exploratory script for local modeling. The new tutorial reports are the better starting point for understanding the ISODISTORT examples themselves.

## Discovery workflow

Use `./.venv/bin/python isotropy_project/scripts/discover_distortion_signatures.py <parent.cif> --label <run_name>` to build a first-pass structural-distortion search from a parent CIF.

Add `--include-coupled-catalog` to also enumerate coupled-primary fixed-k subgroup intersections for the screened vector-active irreps. Use `--max-coupled-pairs` to cap the pair count when the search space becomes too large.

The discovery script now treats the local `iso` output in a way that matches the manuals more closely:

- it enumerates fixed special k-points, irreps, and single-irrep isotropy-subgroup branches locally
- it treats `DISPLAY DISTORTION` rows as projected-vector contributions in the parent display basis
- it transforms those points and vectors into subgroup fractional coordinates before matching them to a constructed child cell
- it collapses repeated contributions on the same child site before applying the probe distortion
- it calls a candidate `verified` only if both the predicted subgroup is recovered and the subgroup-cell embedding is geometrically exact or good
- it can optionally catalog coupled-primary subgroup intersections, which is necessary when no single primary irrep can account for the observed child symmetry
- it records per-query `iso` timing and cache-hit summaries so long runs can be profiled explicitly
- it supports persistent ISO caches plus resumable checkpoints for the parent, k-point, branch, and irrep-screen stages
- it can optionally sample parameterized commensurate k manifolds with explicit `KVALUE` inputs for catalog purposes
- it can optionally catalog one-arm and kernel subgroup counts for additional symmetry context

Each run writes HTML/JSON/Markdown summaries plus candidate CIFs and diffraction plots under `outputs/discovery_runs/<run_name>/`.

The generated reports now include a transition-class coverage table so the workflow states explicitly which cases are implemented, which are only cataloged, and which still require other tools or a different script.

See `docs/isotropy_transition_model.md` for the current conceptual model behind the automation.

For a broader reverse-research summary, environment status, and a proposed roadmap toward a systematic phase-transition explainer, see `docs/systematic_phase_transition_research.md`.

The transition router now also inspects CIF and distortion metadata for magnetic-space-group tags, moment loops, and superspace markers so magnetic and incommensurate cases are dispatched to dedicated workflow families instead of the structural discovery or generic tutorial-report routes.

For the concrete implementation upgrades added after that research pass, see `docs/phase_search_upgrade_notes.md`.

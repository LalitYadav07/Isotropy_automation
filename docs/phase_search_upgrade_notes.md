# Phase Search Upgrade Notes

This note records the concrete implementation upgrades added after the initial reverse-research pass.

## Implemented now

### 1. Discovery instrumentation and auditability

`scripts/discover_distortion_signatures.py` now records:

- per-query `iso` runtime totals
- query call counts
- query cache-hit counts
- total cache-hit count at the run level

These data are written into the JSON report and surfaced in the HTML / Markdown reports so long searches can be profiled explicitly rather than treated as opaque black boxes.

The discovery workflow also now supports persistent ISO query caching and resumable checkpoints so expensive parent standardization, k-point enumeration, branch catalogs, and irrep screens can be reused across reruns.

### 2. `findsym`-aware parent standardization fallback

The new helper `scripts/findsym_tools.py` attempts standardization through the vendored `findsym` workflow and falls back to `pymatgen` conventional standardization only when needed. This gives a more science-aligned route for the parent normalization stage while remaining robust in environments where `findsym` CIF round-tripping is imperfect.

### 3. Transition-family classification

The new module `scripts/phase_transition_classifier.py` classifies the problem before choosing a workflow. The current families are:

- `symmetry_connected_parent_only`
- `symmetry_connected_parent_child`
- `reconstructive_parent_child`
- `decomposition_from_distortion_file`
- `child_only_unsupported`
- `insufficient_input`

The classification is intentionally conservative and reports rationale plus a recommended workflow.

### 4. Top-level workflow router

The new script `scripts/explain_phase_transition.py` acts as a front-end dispatcher:

- parent only -> discovery workflow
- parent + child reconstructive candidate -> `comsubs` wrapper
- distortion file -> tutorial/decomposition report workflow

A `--dry-run` option is included so classification can be validated without launching an expensive downstream search.

### 5. Reconstructive / common-subgroup wrapper

The new script `scripts/reconstructive_transition_search.py` provides an initial parent+child wrapper around `comsubs`. It:

- standardizes parent and child
- prepares a `comsubs` input deck
- runs the vendored executable
- preserves raw output artifacts
- writes a structured JSON summary
- writes a Markdown summary alongside the raw output artifacts

This is deliberately conservative: the parser keeps the raw `comsubs` output as the authoritative artifact and only extracts a minimal machine summary.

### 6. Parameterized-k and one-arm/kernel catalog support

The discovery script now has optional catalog support for:

- sampled parameterized commensurate k manifolds through explicit `KVALUE` inputs
- kernel subgroup counts
- one-arm subgroup counts

These are catalog features only for now. They are not yet carried through to full child-structure realization.

## Scientific cautions that still apply

- Sampled parameterized-k catalogs are not exhaustive. Rational `KVALUE` choices are a controlled probe of the search space, not a proof that all commensurate possibilities have been covered.
- Kernel and one-arm catalogs are symmetry summaries, not verified child structures.
- The reconstructive wrapper is a workflow scaffold around `comsubs`, not yet a complete scientific interpretation layer.
- Single-irrep commensurate vector-displacement discovery remains the most mature part of the local automation.

## Recommended next coding step

The most valuable next engineering step is still to separate discovery into auditable sub-stages with resumable artifacts so that expensive subgroup enumeration does not need to be repeated every time a downstream modeling detail changes.

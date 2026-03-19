# Reverse-Research Notes: Toward a Systematic ISOTROPY Phase-Transition Explainer

This note summarizes what the current repository already does, what was validated in the execution environment, what the ISOTROPY manuals and tutorials imply for a robust workflow, and what should be built next if the goal is a tool that can **systematically understand and explain a structural phase transition**.

## 1. Environment and execution status

### What is installed and working

The Python runtime in this environment already has the three project dependencies available:

- `matplotlib`
- `numpy`
- `pymatgen`

The vendored BYU binaries also run locally once the repository helper configures the process environment. In particular, `scripts/_isotropy_paths.py` creates output directories, sets `ISODATA`, and prepends `vendor/isobyu/` to `PATH`, which is sufficient for `iso`, `findsym`, and related executables to launch.

### What was validated successfully

The following pieces were validated directly:

1. `iso` launches and responds normally.
2. `scripts/generate_tutorial_reports.py` runs successfully and regenerates the tutorial report outputs.
3. Basic k-point enumeration through the lightweight helper scripts works.
4. The discovery workflow successfully parses and standardizes a parent CIF, then enumerates k-points.

### What is not yet fully validated end-to-end

The full `scripts/discover_distortion_signatures.py` workflow is **not yet proven end-to-end in this environment** for the NaCl smoke test.

Observed behavior:

- startup succeeds
- CIF parsing succeeds
- parent standardization succeeds
- `iso` k-point enumeration succeeds
- the full branch enumeration / screening stage takes longer than a 120 s watchdog window on the tested run

This does **not** mean the workflow is conceptually wrong; it means the present environment validation should be described as **partially confirmed, not yet fully performance-qualified**.

## 2. What the repository already understands correctly

The current repo is already built around the correct central idea: the search object is a **symmetry hypothesis**, not an arbitrary set of atom displacements.

The local conceptual model already identifies the key tuple:

- parent space group
- k-star
- irrep
- order-parameter direction
- isotropy subgroup
- subgroup basis and origin

That is the right backbone for a phase-transition explainer because it gives a principled route from the parent crystal to candidate child symmetries.

The repo also already states several critical facts that should be preserved:

1. `vector_active` means only that symmetry allows a microscopic displacement channel on occupied parent sites.
2. `vector_active` does not imply a soft mode or experimental relevance.
3. `DISPLAY DISTORTION` rows are not always one-row-per-atom; they may be projected-vector contributions in the display basis.
4. Those contributions must be transformed into the subgroup basis and merged on the same child-cell site before building a probe structure.
5. Reconstructive transitions are not parent-only discovery problems.
6. Incommensurate, scalar/occupational, and one-arm/kernel workflows remain outside the current automation scope.

These are not minor implementation details. They are the difference between a symmetry-faithful workflow and one that produces attractive but physically misleading subgroup stories.

## 3. What the main scripts do today

### `scripts/discover_distortion_signatures.py`

The current discovery pipeline is already a strong foundation for a future explainer:

1. standardize the parent CIF
2. enumerate fixed special k-points
3. enumerate irreps for each fixed k-point
4. enumerate isotropy subgroups for each irrep
5. retain subgroup basis, origin, and OPD information
6. screen for microscopic vector channels on occupied parent Wyckoff sites
7. build a canonical child cell from subgroup basis and origin
8. transform and aggregate projected distortion rows into child-cell coordinates
9. match the transformed rows onto the child structure
10. apply a small canonical probe distortion
11. independently verify the resulting structure symmetry
12. produce reports and optional coupled-primary catalogs

This is already close to a systematic *single-irrep commensurate displacive* explorer.

### `scripts/generate_tutorial_reports.py`

This script is not a search engine, but it is important for the future product direction. It already behaves more like an **explanatory layer** than a raw computation wrapper:

- it annotates tutorial examples with human-readable intent
- it distinguishes tutorial intent from exported facts
- it reports when ordinary CIF parsing is insufficient for magnetic or incommensurate cases
- it builds HTML / JSON / Markdown outputs appropriate for downstream presentation

In other words: the repo already contains both a **search side** and an **explanation side**. The desired future tool should probably combine them rather than replace either one.

## 4. What the manuals imply for a robust pipeline

The vendored manual and tutorial support several design conclusions.

### 4.1 Coupled order parameters are first-class, not edge cases

The ISOTROPY manual explicitly distinguishes `DISPLAY ISOTROPY` from `DISPLAY ISOTROPY COUPLED`. The coupled display is for subgroups arising from **two or more irreps together**, and the order-parameter direction is shown as a combination of domain-labelled directions from the uncoupled branches.

That means a real transition explainer cannot stop at single-irrep screening. For many experimentally relevant transitions, the child symmetry is the **intersection** of multiple primary order parameters.

### 4.2 One-arm and kernel searches matter strategically

The tutorial notes that long subgroup lists can be shortened by restricting to **one arm of the star** using `VALUE DIRECTION ONEARM`, and that the **kernel** can be queried using `VALUE DIRECTION KERNEL`.

This matters for two reasons:

- performance / search-space control
- interpretability of generic versus symmetry-special OPDs

A systematic explainer should explicitly distinguish:

- full-star isotropy branches
- one-arm branches
- kernel / generic-direction branches

These are not cosmetic variants; they answer different scientific questions.

### 4.3 Continuous-transition information is available and should be surfaced

The manual states that `SHOW CONTINUOUS` indicates which transitions are symmetry-allowed to be continuous in Landau or RG theory.

That is exactly the kind of information users expect from an explainer. Even if the tool does not attempt energetic ranking, it should say things like:

- this subgroup is symmetry-allowed from irrep X
- this OPD is Landau-allowed continuous / not continuous
- this branch is generic versus special

### 4.4 Domains are explanatory output, not optional decoration

`SHOW DOMAINS` gives the domain multiplicity and generators. For a phase-transition explainer, domain structure is central because it connects the abstract subgroup result to experimentally visible variants, twin laws, and switching pathways.

### 4.5 Reconstructive transitions require a different pipeline entirely

The `comsubs` documentation makes clear that it expects **two known crystals** and searches common subgroups under constraints such as cell size, strain, neighbor distance, and shuffle.

This is not just an extension of parent-only irrep enumeration. It is a separate workflow family. Any future tool should therefore begin by classifying the user problem as one of:

- parent-only symmetry-connected discovery
- parent + child symmetry-mode decomposition
- parent + child reconstructive / common-subgroup search

## 5. External research takeaways

Official BYU sources reinforce the same direction:

- the ISOTROPY tutorial emphasizes order-parameter directions, subgroup basis/origin, and domain-aware subgroup generation
- the ISOTROPY introduction paper emphasizes that coupled order parameters can form the primary order parameter of the transition
- ISODISTORT-oriented material emphasizes decomposing structures into symmetry modes and identifying primary versus secondary distortions in a constrained subgroup setting

That triangulates to a clear product goal:

> The right tool is not merely a subgroup enumerator and not merely a mode-decomposition viewer. It is a **decision pipeline** that classifies the transition type, enumerates the correct symmetry hypotheses, and then explains how the observed or candidate child structure follows from primary and secondary order parameters.

## 6. Reverse-engineered product definition

The future tool should answer the following questions for a user:

1. **What kind of transition problem is this?**
   - parent-only discovery
   - parent + child decomposition
   - reconstructive comparison
   - incommensurate superspace problem
   - occupational / scalar ordering problem

2. **What are the candidate primary order parameters?**
   - k-star
   - irrep
   - dimensionality
   - OPD label and vector form
   - generic versus special direction
   - Landau-continuous status when available

3. **What child symmetries follow from each candidate?**
   - subgroup number and symbol
   - basis and origin
   - domain count
   - distinct domain sets where useful

4. **What microscopic content does each candidate imply?**
   - displacive channels on occupied sites
   - strain channels
   - secondary irreps allowed by symmetry
   - whether the branch is only a catalog entry or has a realizable child-cell probe built

5. **How does this compare with a known child structure, if one exists?**
   - exact symmetry match
   - nearly matching subgroup embedding
   - need for coupled-primary explanation
   - need for reconstructive / `comsubs` explanation

6. **How should this be explained in plain language?**
   - what symmetry is broken
   - what remains preserved
   - what the primary OP is doing physically
   - what the secondary modes mean
   - what is established by symmetry and what is not

## 7. Proposed systematic pipeline

## Stage A — problem classification

Input options:

- parent CIF only
- parent CIF + child CIF
- parent CIF + ISODISTORT distortion file
- parent CIF + experimental constraints (known child SG, known modulation vector, known active atom set, etc.)

Decision tree:

1. if a distortion file is provided, prefer the **explanation / decomposition** path
2. if both parent and child are provided, decide whether the phases are symmetry-connected or reconstructive
3. if only a parent CIF is provided, use the **discovery** path
4. if the k vector is incommensurate or parameterized, route away from the current commensurate builder and into a dedicated superspace branch

## Stage B — parent normalization

Use a more robust standardization layer than the current `pymatgen`-only approach when possible.

Target behavior:

- conventional standardized parent cell
- verified parent SG and setting
- occupied Wyckoff list
- explicit mapping from input cell to standardized cell

A direct scripted `findsym` pathway should likely become the canonical standardization route, with `pymatgen` retained as a fallback or convenience layer.

## Stage C — symmetry hypothesis generation

For parent-only discovery, enumerate in this order:

1. fixed special k-points
2. irreps at each fixed k-point
3. single-irrep isotropy subgroups
4. one-arm variants when relevant
5. kernel / generic-direction variants when relevant
6. coupled-primary subgroup intersections among screened irreps
7. parameterized commensurate k-lines / k-planes under explicit rational sampling policies

This stage should emit a **deduplicated hypothesis table**, not just raw command output.

## Stage D — structural realizability filters

For each hypothesis, compute or annotate:

- occupied-site vector activity
- occupied-site scalar / occupational activity
- child-cell size and basis determinant
- whether a canonical child-cell probe can be built
- whether a known child structure matches the subgroup embedding
- whether the case is catalog-only because full structure generation is not implemented

This is the stage where “symmetry allowed” is separated from “structurally instantiated in this automation.”

## Stage E — child comparison / decomposition

If a child structure exists:

1. compare parent and child after standardization
2. test whether the child can be explained by a single-irrep branch
3. if not, test coupled-primary subgroup intersections
4. if still not explained, route to reconstructive / `comsubs`
5. if an ISODISTORT distortion export exists, parse it and summarize primary versus secondary channels directly

## Stage F — explanation synthesis

For each plausible pathway, generate:

- concise scientific summary
- symmetry chain in plain language
- order-parameter explanation
- active secondary modes
- domain count and meaning
- why the branch is plausible or implausible given the provided evidence
- what remains unknown without energetics, phonons, or experiment

This is the stage that converts a technically correct catalog into a real explainer.

## 8. Recommended near-term implementation roadmap

### Priority 1 — make discovery performance auditable

Before expanding scope, make the current discovery runtime transparent.

Add instrumentation for:

- time spent in each `iso` call
- counts of irreps and subgroup branches per k-point
- cache hits / misses
- branch-screening throughput
- reasons for long-running steps

Without this, it is difficult to know whether the current slowness is due to expected combinatorics, repeated identical queries, or command-level inefficiency.

### Priority 2 — add a first-class “transition classifier” layer

Implement a small front-end module that classifies the problem before choosing the computational path.

Suggested outputs:

- `symmetry_connected_parent_only`
- `symmetry_connected_parent_child`
- `reconstructive_parent_child`
- `incommensurate_superspace`
- `occupational_or_scalar`
- `decomposition_from_distortion_file`

### Priority 3 — promote coupled-primary analysis from catalog to explanation

The repo already catalogs coupled-primary fixed-k subgroup intersections. The next step is to make coupled-primary results explainable:

- build coupled child structures where feasible
- verify them like the single-irrep probes
- explain which irreps act jointly as the primary OP
- report when the child symmetry cannot be obtained from any single primary irrep but is recovered as an intersection

### Priority 4 — add parent + child reconstructive workflow around `comsubs`

This is the clearest missing workflow family and aligns with the existing project notes.

The wrapper should:

- standardize both parent and child
- extract required `comsubs` inputs
- expose sensible constraints (`size`, `strain`, `neighbor`, `shuffle`, optional charges)
- summarize candidate common-subgroup pathways in a human-readable report

### Priority 5 — expand beyond vector-only screening

A structural transition explainer should eventually cover:

- scalar / occupational order
- rotational order where distinct from simple displacive summaries
- incommensurate superspace cases
- parameterized commensurate manifolds with rational sampling

## 9. Suggested architecture for the future tool

A clean architecture would likely have five layers:

1. **Normalization layer**
   - CIF parsing
   - standardization
   - parent/child mapping

2. **Symmetry query layer**
   - `iso`
   - `findsym`
   - `comsubs`
   - caching and timing

3. **Structural realization layer**
   - child-cell construction
   - distortion transformation / aggregation
   - verification and comparison

4. **Transition reasoning layer**
   - classify transition family
   - rank explanatory pathways
   - identify single-irrep, coupled, or reconstructive interpretations

5. **Narrative/report layer**
   - JSON for machine use
   - Markdown / HTML for people
   - concise explanatory summaries
   - auditable evidence sections

## 10. Bottom line

The project is already pointed in the right direction.

The current codebase is **not** just a random set of wrappers; it already contains the key intellectual move needed for a serious ISOTROPY-based transition analyzer: treat the transition as a search over symmetry hypotheses defined by irreps, order-parameter directions, and isotropy subgroups.

The main gap is not the absence of a good idea. The main gaps are:

- incomplete workflow coverage beyond fixed-k single-irrep commensurate displacive cases
- lack of a top-level transition classifier
- coupled-primary results are cataloged but not yet fully realized and explained
- reconstructive analysis is recognized but not yet wrapped
- incommensurate and scalar/occupational cases are not automated
- discovery runtime needs better instrumentation

If those gaps are addressed in order, this repository can evolve into a tool that does what you want: **systematically understand and explain structural phase transitions using ISOTROPY rather than merely listing subgroups.**

## 11. Useful official references consulted

- ISOTROPY tutorial PDF: <https://iso.byu.edu/isotut.pdf>
- ISOTROPY introduction / overview PDF: <https://iso.byu.edu/iso/2006%20Stokes.pdf>
- Local vendored manuals used during reverse research:
  - `vendor/isobyu/isoman.txt`
  - `vendor/isobyu/isotut.txt`
  - `vendor/isobyu/comsubs.txt`

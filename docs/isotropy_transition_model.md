# ISOTROPY Transition Model

This note records the symmetry concepts that the local automation is trying to respect.

## The right search object

For symmetry-connected structural transitions, the primary search object is not an arbitrary list of atomic shifts. It is a symmetry hypothesis:

- parent space group
- k-star
- irrep
- order-parameter direction
- isotropy subgroup
- subgroup basis and origin

For coupled transitions, the search object becomes an intersection of two or more uncoupled isotropy subgroups. In ISOTROPY language this is a coupled order parameter, not a single-irrep branch.

## What single-irrep isotropy search means

If a single primary irrep is sufficient, the isotropy subgroup already determines:

- the child symmetry
- the allowed secondary order parameters
- the allowed strain couplings
- the domain structure

That is why the first brute-force layer should be a catalog of unique subgroup hypotheses, not arbitrary mode-amplitude combinations.

## Structural transition classes

The main classes relevant to this project are:

1. Commensurate fixed-k displacive transitions
   These are the cleanest cases for local automation with `iso` plus a canonical probe distortion.

2. Commensurate fixed-k coupled-primary transitions
   These appear when no single irrep explains the observed child symmetry. They must be treated as coupled order parameters or subgroup intersections.

3. Parameterized commensurate k-line and k-plane transitions
   These involve k vectors with free parameters and need an additional search over commensurate choices or experimental constraints.

4. Incommensurate transitions
   These require superspace treatment. A conventional 3D child-cell builder is not enough.

5. Occupational or order-disorder transitions
   These can involve microscopic scalar distortions rather than atomic displacement vectors.

6. Reconstructive transitions
   These are not captured by a single group-subgroup chain from the parent. They require parent-child comparison and common-subgroup analysis, for example with `comsubs`.

## Practical interpretation

- `findsym` standardizes and verifies structures.
- `iso` enumerates irreps, isotropy subgroups, domains, and distortions for symmetry-connected branches.
- `DISPLAY ISOTROPY COUPLED` is needed when the child symmetry requires coupled primary order parameters.
- `comsubs` is needed when both phases are known but the transition is reconstructive.

## Current local workflow

The current discovery workflow in `scripts/discover_distortion_signatures.py` does the following:

- standardizes the parent CIF
- catalogs fixed-k single-irrep subgroup branches
- screens fixed-k irreps for vector-active structural channels on occupied Wyckoff sites
- models primary single-irrep vector distortions as small canonical probes
- catalogs coupled-primary subgroup intersections for screened vector-active irrep pairs when requested
- reports explicitly which transition classes are implemented, cataloged only, or still out of scope

Here "vector-active" means only that symmetry allows a microscopic displacement pattern on the occupied parent Wyckoff sites. It does not mean the irrep is dynamically unstable, energetically favorable, or experimentally present.

## Common failure mode to avoid

`DISPLAY DISTORTION` does not necessarily mean one row equals one distinct atom. Depending on the case, rows can represent projected-vector contributions in the chosen display cell. Those contributions must be interpreted in the correct basis and combined on the same child-cell site before using them to build a trial structure.

That point is essential. If it is mishandled, the workflow can appear to "verify" or "reject" subgroups for purely geometric reasons that have nothing to do with the underlying symmetry physics.


## Known magnetic child structures

When both a parent CIF and a known magnetic child CIF are available, the safest automated path is now a decomposition-centered workflow rather than parent-only discovery. In practice that means: standardize both settings, keep the saved ISODISTORT basis/origin/OPD metadata explicit, treat the dominant magnetic mode family as the primary order-parameter candidate, and report non-magnetic amplitudes as secondary displacive / ordering / rotational / strain responses within the same subgroup scaffold.

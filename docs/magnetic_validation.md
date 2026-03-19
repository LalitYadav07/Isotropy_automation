# Magnetic tutorial validation fixtures

This document turns the magnetic tutorial artifacts under `isotutorials/exercises-completed/` into explicit regression fixtures. The goal is to prevent future parser or reporting changes from silently breaking magnetic or superspace support.

## Scope and intent

These fixtures cover the magnetic tutorial examples already used by `scripts/generate_tutorial_reports.py`:

- `lamno3-magchild-distortion.txt`
- `dymn6ge6-distortion.txt`
- `tbmno3-distortion.txt`
- `tbmno3-distortion2.txt`
- `skyrmion-distortion.txt`

Each fixture records:

1. the parent structure artifact,
2. the active magnetic irreps taught by the tutorial,
3. the OPD / subgroup selection to preserve,
4. expected primary and secondary channels, and
5. the exact facts local code should extract or verify.

The machine-readable source of truth lives in `fixtures/magnetic_tutorial_regression.json`, and `tests/test_magnetic_tutorial_regression.py` enforces it.

## Regression matrix

| Example | Parent structure | Active magnetic irreps | OPD / subgroup | Expected primary channels | Expected secondary channels | What local code must extract or verify |
| --- | --- | --- | --- | --- | --- | --- |
| LaMnO3 mode decomposition | `lamno3-cubic.cif` | `mX5+` tutorial primary order | `Subgroup: 62.448 Pn'ma'` with basis `{(1,0,-1),(0,2,0),(1,0,1)}` | Antiferromagnetic `mX5+` order | Nonzero displacive and strain decomposition terms | Recover the stored magnetic subgroup, confirm magnetic/displacive/strain families, and confirm that this decomposition export does **not** expose explicit primary-irrep labels in the parsed `irrepString` field. |
| DyMn6Ge6 conical order | `dymn6ge6-parent.cif` | `mGM2+`, `mDT6` | `P1-P (a|b,0,0,0)` → `177.1.24.2.m153.1 P62'2'(0,0,g)h00` | Ferromagnetic `mGM2+` plus incommensurate `mDT6` cone | Strain scaffold present but zero | Extract both primary irreps, both k vectors, all seven nonzero magnetic coefficients, and preserve the superspace subgroup string. |
| TbMnO3 cycloid, primary only | `tbmno3-pbnm.cif` | `mSM3`, `mSM2` | `P-P (a,0|b,0)` → `33.1.9.5.m145.2 Pbn2_1.1'(0,a,0)000s` | Cycloidal magnetic pair `mSM3 + mSM2` | Polar displacements symmetry-allowed but zero | Extract the two magnetic irreps, confirm only magnetic amplitudes are populated, and keep the polar/displacive scaffold visible. |
| TbMnO3 cycloid with polar response | `tbmno3-pbnm.cif` | `mSM3`, `mSM2` | same subgroup as the primary-only file | Cycloidal magnetic pair `mSM3 + mSM2` | Explicit secondary polar displacements | Extract the same magnetic subgroup and irreps as above, but also detect the two nonzero displacive coefficients. |
| Fe-monolayer skyrmion scaffold | `skyrmion-parent.cif` | `mLD3`, `mLD4` | `P-P (a,0;a,0;a,0|b,0;b,0;b,0)` → `177.2.83.6.m153.1 P62'2'(a,a,0)000(-2a,a,0)000` | Multi-k magnetic skyrmion/vortex scaffold | Occupational, displacive, and strain channels allowed but zero | Preserve the two magnetic irreps, the three-arm `k-active` list, and the zero-valued secondary-channel scaffold instead of dropping it during parsing. |

## Design notes

### 1. Tutorial truth vs parser-observable truth

The fixtures intentionally separate **tutorial interpretation** from **parser observables**.

- For `lamno3-magchild-distortion.txt`, the tutorial teaches `mX5+` as the primary magnetic order.
- The exported decomposition file, however, stores only the chosen subgroup plus mode coefficients and does not expose the primary irrep labels through `irrepString`.

The regression fixture therefore records both facts: the tutorial-level magnetic irrep and the parser-level expectation that `primary_irreps == []`.

### 2. Zero-valued scaffolds are still important

Several magnetic examples are useful precisely because they preserve symmetry scaffolding even when coefficients are zero.

- `tbmno3-distortion.txt` keeps a polar/displacive scaffold with zero amplitudes.
- `skyrmion-distortion.txt` keeps magnetic, occupational, displacive, and strain metadata even though the stored coefficients are zero.

Future code should not prune those zero-valued sections away if the goal is faithful magnetic support.

### 3. What the regression test actually guards

The automated test asserts that local parsing continues to recover:

- distortion classification labels,
- exact or substring-stable OPD strings,
- expected primary irrep and k-vector records where the file exposes them,
- nonzero mode counts by family,
- child-structure atom counts, and
- parser internal consistency checks.

That means future changes which accidentally drop magnetic metadata, scramble irrep ordering, or erase superspace OPD text will fail fast.

## How to run the regression check

```bash
python -m unittest tests.test_magnetic_tutorial_regression
```

If we later add dedicated magnetic extraction logic beyond `parse_distortion_file`, this fixture set should become the baseline acceptance suite for that code as well.

# Glycine-specific Voronoi collective variables

This directory archives the PLUMED C++ collective variables used for the glycine tautomerism workflows in this repository. The implementation is differentiable and supports MPI/OpenMP execution, but it is a **paper-specific legacy interface**, not the current general Reactive Soft-Voronoi API.

For the general design, explicit atom-selection rules, current API, and the recommended validation workflow, see the [Reactive Soft-Voronoi collective-variable guide](https://zhang-pchao.github.io/code/reactive-voronoi). The guide explains the common smooth assignment and its use for glycine tautomerism, bulk-water autoionization, and interfacial ion location.

## Contents

| Source | Registered PLUMED Action | Role in the archived workflow |
| --- | --- | --- |
| [`VoronoiC0.cpp`](./VoronoiC0.cpp) | `VORONOIC0` | Squared coordination-defect activity over the water-like sites; used as the ion/charge-activity coordinate. |
| [`VoronoiD2.cpp`](./VoronoiD2.cpp) | `VORONOID2` | Defect-weighted distance between water sites and the glycine reactive sites; used as the defect-separation coordinate. |
| [`VoronoiS1.cpp`](./VoronoiS1.cpp) | `VORONOIS1` | Weighted protonation/solvation coordinate combining the water and glycine site defects. |

The repository-level examples that define the published atom groups are [`Enhanced_Sampling/config_sample/plumed.dat`](../Enhanced_Sampling/config_sample/plumed.dat) and [`research/GlycineEfield/Enhanced_Sampling_MD/LongRange_Interface/input.plumed`](https://github.com/Zhang-pchao/research/blob/main/GlycineEfield/Enhanced_Sampling_MD/LongRange_Interface/input.plumed).

## Legacy data model

All three Actions use the same underlying pattern:

1. `GROUPA` contains candidate centers.
2. `GROUPB` contains transferable atoms, such as water H, N--H, and carboxyl O--H atoms.
3. Each `GROUPB` atom is smoothly assigned to the `GROUPA` centers by a distance-weighted kernel.
4. `D_0`, `D_1`, `D_2`, and `D_3` define reference occupancies by **position in `GROUPA`**.
5. The final Action reduces the resulting coordination defects into a scalar CV.

For the glycine examples, `GROUPA` must be ordered as:

```text
water O centers, then glycine N, glycine O1, glycine O2
```

and `NRX=3` tells the code that the final three entries are the special glycine sites. This positional convention is part of the legacy implementation. Reordering `GROUPA` without simultaneously changing `NRX` and `D_0...D_3` changes the CV definition.

The source evaluates the kernel as `exp(LAMBDA * distance)`. The archived inputs therefore use negative values such as `LAMBDA=-5`, `-8`, or `-100`. Do not interpret `LAMBDA` as the positive `KAPPA` parameter in the current guide; scan and validate the value on representative structures before using it for biasing.

## Build a runtime plugin

Use the same PLUMED executable, compiler, MPI mode, and ABI as the target MD engine. A convenient development path is PLUMED's runtime-library builder:

```bash
cd /path/to/GlycineTautomerism/Voronoi_collective_variables
mkdir -p build-voronoi-glycine
cd build-voronoi-glycine

plumed mklib ../VoronoiC0.cpp
plumed mklib ../VoronoiD2.cpp
plumed mklib ../VoronoiS1.cpp
```

If the local PLUMED installation supports compiling several source files in one `mklib` invocation, they may be built as one plugin instead. The important contract is that the resulting `.so` files are compiled against the same PLUMED installation used by the simulation. Rebuild after changing PLUMED, the compiler, MPI, or the ABI.

Load the libraries before the Actions are referenced:

```plumed
LOAD FILE=./VoronoiC0.so
LOAD FILE=./VoronoiD2.so
LOAD FILE=./VoronoiS1.so
```

Alternatively, copy the sources into the matching `plumed/src/colvar` tree and rebuild PLUMED. Runtime loading is preferable for an initial audit because it keeps the experiment isolated from the main PLUMED installation.

## Minimal archived input

The following is a template, not a universal parameter set. Replace the groups with the atom numbering in the coordinate file read by PLUMED.

```plumed
LOAD FILE=./VoronoiC0.so
LOAD FILE=./VoronoiD2.so
LOAD FILE=./VoronoiS1.so
UNITS LENGTH=A

# GROUPA order: water O, glycine N, glycine O1, glycine O2
WaterO:    GROUP ATOMS=...
GlyN:      GROUP ATOMS=...
GlyO1:     GROUP ATOMS=...
GlyO2:     GROUP ATOMS=...
WaterH:    GROUP ATOMS=...
GlyH:      GROUP ATOMS=...
Centers:   GROUP ATOMS=WaterO,GlyN,GlyO1,GlyO2
AllH:      GROUP ATOMS=WaterH,GlyH

# NRX=3: the last three GROUPA entries are GlyN, GlyO1, and GlyO2.
sd: VORONOID2 GROUPA=Centers GROUPB=AllH NRX=3 LAMBDA=-5 \
    D_0=2 D_1=2 D_2=0.5 D_3=0.5
sp: VORONOIS1 GROUPA=Centers GROUPB=AllH NRX=3 LAMBDA=-5 \
    D_0=108 D_1=3
sa: VORONOIC0 GROUPA=Centers GROUPB=AllH NRX=3 LAMBDA=-5 \
    D_0=2 D_1=2 D_2=0.5 D_3=0.5

PRINT ARG=sd,sp,sa FILE=COLVAR STRIDE=1
DUMPDERIVATIVES ARG=sd,sp,sa FILE=DERIVATIVES STRIDE=1
```

`VORONOIS1` registers only `D_0` and `D_1`; do not pass `D_2` or `D_3` to that Action. The `D_0=108`, `D_1=3` values above match the archived 54-water glycine example and must be recomputed when the number of centers or the reference chemical state changes.

## Neighbor-list and validation contract

With no `NLIST`, the code evaluates the full available center--assigned pair set and this should be the reference calculation. When `NLIST` is enabled, both `NL_CUTOFF` and `NL_STRIDE` are required:

```plumed
sd_fast: VORONOID2 GROUPA=Centers GROUPB=AllH NRX=3 LAMBDA=-5 \
    D_0=2 D_1=2 D_2=0.5 D_3=0.5 \
    NLIST NL_CUTOFF=2.4 NL_STRIDE=1
```

The `2.4` value is an archived water-system setting, not a default or transferable bond cutoff. This legacy code does not implement the newer `NL_SKIN` keyword. Use the following sequence before production sampling:

1. Run exact mode on labeled neutral, transition, product, host-switching, and distorted frames.
2. Compare the accelerated value **and coordinate/box derivatives** with exact mode while increasing `NL_CUTOFF`.
3. Test `NL_STRIDE>1` only after the cutoff has converged; for replica exchange, choose a stride compatible with the exchange schedule.
4. Use `DUMPDERIVATIVES` and a numerical-derivative diagnostic on a short fixture. A successful `plumed driver` run checks parsing and evaluation, not the correctness of a bias or production force.
5. Record the PLUMED version, source commit, compiler/MPI ABI, `LAMBDA`, reference values, NLIST settings, comparison tolerances, and checksums of the plugin files.

## Current API migration

For new systems, the current guide is the preferred engineering interface:

| Archived Action | Current conceptual replacement |
| --- | --- |
| `VORONOIC0` | `VORONOI_COORDINATION ... POWER=2` |
| `VORONOIS1` | `VORONOI_COORDINATION ... POWER=1`, with explicit coefficients and selections |
| `VORONOID2` | `VORONOI_DISTANCE` with explicit `GROUP1/GROUP2` pair semantics |

The current API removes the positional `NRX` convention and declares `CENTERS`, `ASSIGNED`, and `REFERENCE` explicitly. Keep this directory for reproducing the paper-specific workflow; use the [full guide](https://zhang-pchao.github.io/code/reactive-voronoi) and its current source/regression examples when designing a new CV.

## Related references

- [Glycine tautomerism repository](https://github.com/Zhang-pchao/GlycineTautomerism)
- [Published glycine study](https://doi.org/10.1021/acs.jcim.4c00273)
- [Archived production input](../Enhanced_Sampling/final_MD/input.plumed)
- [Current Reactive Soft-Voronoi guide](https://zhang-pchao.github.io/code/reactive-voronoi)

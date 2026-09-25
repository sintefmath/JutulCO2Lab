# JutulCO2Lab

Exploratory reimplementation of MRST-co2lab's vertical-equilibrium (VE) CO2
storage simulator on top of [Jutul](https://github.com/sintefmath/Jutul.jl) /
[JutulDarcy](https://github.com/sintefmath/JutulDarcy.jl).

What exists: the environment, the oracle harness, the VE domain (Stage 1:
columns, geometry, rock, caprock datum, boundary sides — all derived from the
3D grid and verified against MRST), and Stage 2 items 0–6: co2lab's fluid
system, VE relperm, `pcWG`, PVT, VE wells, VE boundary conditions and
reconstruction of 3D states, running end to end on co2lab's sloping aquifer.
A readable walkthrough is `../scripts/example_sloping_aquifer.jl`. The plan
lives in `../analysis/roadmap.org`; the evidence behind it in
`../analysis/ve-meshing/hybrid-ve-feasibility.md`.

## Environment

`Jutul` and `JutulDarcy` are `dev`-ed from sibling checkouts, not from the
registry — the port tracks those working trees, and several roadmap decisions
cite them by `file:line`:

- `Jutul`: `../../Jutul.jl-fork`, the researcher's fork (github.com/Oddan/Jutul.jl)
  checked out on branch `co2lab-jutul-merge` = the pinned upstream commit
  `36166ba6` plus the Jutul fixes this project needs (so far: the cross-term
  row stride). `../../Jutul.jl` (upstream) is kept at `36166ba6` for reference.
- `JutulDarcy`: `../../JutulDarcy.jl`.

For a dev-ed package, `Manifest.toml` records only `path` and `version` — no
commit — so it does **not** pin what actually ran. Every `OracleReport`
therefore records the commit and dirty state of both checkouts when it is
created, and prints them in its header:

```
Stage 1 — Johansen grid
  Jutul 36166ba6  JutulDarcy 2c92d4d5  MRST 1a2b3c4d
```

Pass the MRST dump as `OracleReport(title; mrst = ref)` to include the MRST
commit. A checkout with uncommitted changes to tracked files triggers a
warning, since the report then cannot be reproduced from its commit alone.
Record-only: nothing fails on provenance. `source_provenance()` gives the same
information on demand.

```
../scripts/julia --project=. -e 'using Pkg; Pkg.instantiate()'
../scripts/julia --project=. -e 'using Pkg; Pkg.test()'
```

The test suite ends with a **toolchain smoke test** that assembles and solves a
small stock immiscible case. If it passes, Jutul/JutulDarcy work here; failures
in it are environment problems, not VE problems.

## The oracle

Every acceptance gate in the roadmap is a numerical comparison of this package
against MRST-co2lab **on identical input** — not a plot, and not a comparison
against a 3D reference. Every trap in feasibility §2b fails silently, so a gate
is a per-quantity check with a declared tolerance.

```julia
using JutulCO2Lab.Oracle

ref = read_mrst_dump("stage1_johansen.mat")
rep = OracleReport("Stage 1 — Johansen grid"; mrst = ref)
check!(rep, "column height", H, ref["height"]; tol = STAGE1_TOL[:height], perm = c2m)
check!(rep, "pore volume",   pv, ref["pore_volume"]; tol = STAGE1_TOL[:pore_volume], perm = c2m)
display(rep)
passed(rep) || error("Stage 1 gate failed")
```

Two deliberate omissions in `Oracle`:

- **It never guesses orderings.** MRST and Jutul cell/face numbering need not
  agree; pass an explicit `perm` (Julia index → MRST index). A silently
  permuted comparison is the exact false pass the harness exists to prevent.
- **An empty report is a failure, not a pass.** A gate that ran no checks has
  not been executed.

The MATLAB half is `../scripts/mrst_dump_stage1.m`, which mirrors
`exampleVEJohansen.m` and `CO2VEBlackOilTypeModel.setupOperators` line for
line. Run it from a MATLAB session with MRST on the path:

```matlab
run /path/to/mrst/startup.m
mrst_dump_stage1('stage1_johansen.mat')
```

It dumps what the Stage 1 gate needs: column heights, pore volumes, the
`T_all` transmissibilities (interior **and** boundary), face connectivity,
`transMult`, the upscaled rock, the ordered `Gt.columns` structure, and the
height-scaled well index.

## Layout

```
src/JutulCO2Lab.jl   module root
src/oracle.jl        Oracle — MRST dump reader and the tolerance-checked comparison
src/provenance.jl    source commits + dirty state, recorded in every report
src/vegrid.jl        VE columns on the 3D mesh: partition, trimming, column geometry
src/vefaces.jl       face directions (MRST tags / logical indices) and column sides — topology only
src/vedomain.jl      VE domain on a CoarseMesh: upscaled rock, caprock datum, boundary sides
src/ve_fluid.jl      co2lab's VE fluid: plume heights, sharp-interface relperm, pcWG,
                     linear PVT and pore volume
src/co2lab_fluid.jl  co2lab's fluid system: black-oil with a :co2lab formulation,
                     unknowns by configuration; model, state and simulator setup
src/vewells.jl       VE wells: co2lab's well index, pressure datum option, CO2 injector
src/vebc.jl          VE boundary conditions: co2lab's pressure BC as a force, side
                     selection, exited mass
src/vereconstruct.jl 3D states from VE states: saturation from plume heights,
                     hydrostatic pressure
test/runtests.jl     Oracle, Stage 1 (vs MRST dumps), fluid system, VE fluid, wells,
                     boundary conditions and reconstruction (vs MRST), toolchain smoke test
../scripts/          MRST/Octave reference dumps, prototypes, gate scripts
```

All types and returned objects, in order of dependency: `STRUCTURES.md`.

## Open decisions

Tracked in `../analysis/STATUS.md`.

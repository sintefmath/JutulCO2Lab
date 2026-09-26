# JutulCO2Lab

Reimplementation of MRST-co2lab's vertical-equilibrium (VE) CO2 storage
simulator on top of [Jutul](https://github.com/sintefmath/Jutul.jl) /
[JutulDarcy](https://github.com/sintefmath/JutulDarcy.jl).

What exists: the environment, the oracle harness, the VE domain (Stage 1:
columns, geometry, rock, caprock datum, boundary sides — all derived from the
3D grid and verified against MRST), and Stage 2: co2lab's fluid system, VE
relperm, `pcWG`, PVT, VE wells, VE boundary conditions and reconstruction of
3D states, running end to end on co2lab's sloping aquifer, plus the
structures Stage 8 (hybrid) extends — stackable columns cut by connectivity
and per-cell VE/fine branch points; Stage 3: co2lab's hysteresis (residual
trapping); and Stage 4: dissolution, rate-driven and instantaneous, each
gated against co2lab to solver precision.
Readable walkthroughs: `examples/sloping_aquifer.jl` (Stage 2),
`examples/residual_trapping.jl` (Stage 3, residual trapping) and
`examples/dissolution.jl` (Stage 4, rate-driven and instantaneous
dissolution). The plan
lives in `../analysis/roadmap.org`; the evidence behind it in
`../analysis/ve-meshing/hybrid-ve-feasibility.md`.

## Architecture

A VE model is a JutulDarcy reservoir model whose grid is a `CoarseMesh` of
VE columns over the 3D grid, and whose black-oil system carries co2lab's
formulation and VE physics as replaced variables. Nothing in Jutul or
JutulDarcy is patched except one Jutul bug fix (see Environment).

**Data flow** (`examples/sloping_aquifer.jl` follows it step by step):

1. *3D mesh*: a Jutul `UnstructuredMesh` with logical `(i, j, k)` and a
   direction tag on every boundary face — from MRST (`ve_mesh_from_mrst`) or
   from Jutul's own corner-point processing.
2. *Columns*: `ve_columns` → `VEColumns`, the partition of fine cells into
   columns (vertically connected stretches of a pillar); `ve_trim` removes
   discarded cells.
3. *VE domain*: `ve_domain` → a `CoarseMesh` over the 3D mesh, wrapped by
   `reservoir_domain`, with thickness-averaged rock, the caprock datum as
   `:gdz`, and `:column_height`, `:rugosity`, `:is_ve`, `:pillar_columns`.
   Transmissibility and pore volume get their ×H from the coarse geometry.
4. *Wells and boundaries*: `setup_ve_well` (co2lab's well index, perforation
   at the caprock); `ve_boundary_sides` → `ve_select_sides` →
   `ve_boundary_conditions` (co2lab's pressure BC, one `VEBoundaryCondition`
   per open side).
5. *Model*: `setup_co2lab_model` → a JutulDarcy `MultiModel` (reservoir, wells,
   facility); `setup_co2lab_state`, `setup_co2lab_simulator`; forces from
   JutulDarcy's `setup_reservoir_forces(...; bc)`; run with JutulDarcy's
   `simulate_reservoir`.
6. *Results*: JutulDarcy's states and well results; `ve_boundary_mass_fluxes`
   (exited mass); `ve_reconstruct` (3D pressure and saturations).

**Where co2lab's physics enters JutulDarcy:**

| Hook | What is replaced or added | File |
|---|---|---|
| system type | `StandardBlackOilSystem` re-tagged with formulation `:co2lab` (`CO2LabSystem`); brine = liquid, CO2 = vapour | `co2lab_fluid.jl` |
| `select_primary_variables!` | unknowns `p`, `sG` (and `rs` with dissolution; `sGmax` with hysteresis, added by `setup_co2lab_model`) | `co2lab_fluid.jl` |
| `select_secondary_variables!` and `@jutul_secondary` | saturations, `TotalMasses`, dissolved CO2; plume heights, sharp-interface relperm, `pcWG`, linear PVT, pore volume | `co2lab_fluid.jl`, `ve_fluid.jl`, `rock_compressibility.jl` |
| parameters (`default_parameter_values`) | `ColumnHeight`, `Rugosity` from the domain; without hysteresis `sGmax`, updated after each converged step (`update_after_step!`) | `ve_fluid.jl`, `co2lab_fluid.jl` |
| extra equations | `DissolvedCO2Equation` (rate-driven dissolution) and `HysteresisEquation` (`sGmax` with hysteresis), each with its own convergence criterion | `co2lab_fluid.jl` |
| `update_primary_variable!` | co2lab's switching unknown for instantaneous dissolution (`DissolvedOrFreeCO2`, on JutulDarcy's `BlackOilX`) | `instant_dissolution.jl` |
| cross term | `DissolutionWellSupplyCT`: a well's CO2 as dissolution supply in its cells | `co2lab_fluid.jl` |
| `apply_forces_to_equation!` | `VEBoundaryCondition` passed as JutulDarcy's `bc` force | `vebc.jl` |
| flux | none: JutulDarcy's stock two-point flux with the caprock `:gdz` is co2lab's VE flux (roadmap Stage 0) | — |
| wells | none: JutulDarcy's wells with co2lab's well index and a CO2 wellbore column (`co2_injector_control`) | `vewells.jl` |

The VE variables branch per cell on `is_ve`, and fluxes will branch per face
through a custom `kgrad` type: these are where Stage 8 (hybrid) adds fine 3D
cells (`../analysis/ve-meshing/isve-branch-proposal.md`).

Deliberate differences from co2lab are listed in the roadmap's standing
decisions ("Deliberate deviations from co2lab").

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

src/ve/              the VE simulator
  vegrid.jl          VE columns on the 3D mesh: partition, trimming, column geometry
  vefaces.jl         face directions (MRST tags / logical indices) and column sides
                     — topology only
  vedomain.jl        VE domain on a CoarseMesh: upscaled rock, caprock datum,
                     boundary sides
  rock_compressibility.jl  co2lab's linear pore volume, a model setting
  ve_fluid.jl        co2lab's VE fluid: plume heights, sharp-interface relperm, pcWG,
                     linear PVT; column height and rugosity parameters
  co2lab_fluid.jl    co2lab's fluid system: black-oil with a :co2lab formulation,
                     unknowns by configuration; dissolution transfer, hysteresis;
                     model, state and simulator setup
  instant_dissolution.jl  co2lab's instantaneous dissolution: switching unknown
  vewells.jl         VE wells: co2lab's well index, pressure datum option, CO2
                     injector
  vebc.jl            VE boundary conditions: co2lab's pressure BC as a force, side
                     selection, exited mass
  vereconstruct.jl   3D states from VE states: saturation from plume heights,
                     hydrostatic pressure

src/atlas/           the CO2 Storage Atlas (Stage 5)
  atlas.jl           data location, map reader, formation names, average rock
  atlas_grid.jl      top and thickness maps to a corner-point grid section, rock maps
  atlas_domain.jl    a formation as a VE domain: mesh, largest piece, columns, rock

examples/            readable walkthroughs (run from the project root with
                     scripts/julia --project=JutulCO2Lab JutulCO2Lab/examples/...)
test/runtests.jl     Oracle, Stage 1 (vs MRST dumps), fluid system, VE fluid, wells,
                     boundary conditions and reconstruction (vs MRST), Atlas,
                     toolchain smoke test
../scripts/          MRST/Octave reference dumps, prototypes, gate scripts
```

All types and returned objects, in order of dependency: `STRUCTURES.md`.

## Open decisions

Tracked in `../analysis/STATUS.md`.

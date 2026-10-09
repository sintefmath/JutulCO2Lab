# JutulCO2Lab

Reimplementation of MRST-co2lab's vertical-equilibrium (VE) CO2 storage
simulator on top of [Jutul](https://github.com/sintefmath/Jutul.jl) /
[JutulDarcy](https://github.com/sintefmath/JutulDarcy.jl).

What exists: the environment, the oracle harness, the VE domain (Stage 1:
columns, geometry, rock, caprock datum, boundary sides — all derived from the
3D grid and verified against MRST), and Stage 2: co2lab's fluid system, VE
relperm, `pcWG`, PVT, VE wells, VE boundary conditions and reconstruction of
3D states, running end to end on co2lab's sloping aquifer, plus the
structures Stage 12 (hybrid) extends — stackable columns cut by connectivity
and per-cell VE/fine branch points; Stage 3: co2lab's hysteresis (residual
trapping); and Stage 4: dissolution, rate-driven and instantaneous, each
gated against co2lab to solver precision; Stage 5: the CO2 Storage Atlas
and the Sleipner benchmark; Stage 8: co2lab-spillpoint's trapping analysis
(both methods, the derived quantities, IRAP/IGEMS surfaces), identical to
co2lab (faults sealing by default); Stage 9 items 1–2: the trapping
inventory, identical to co2lab on identical input and, end to end, against
co2lab's own run of a Johansen case; Stage 6 (in progress): the
vertical-model framework (`VerticalModel`; co2lab's sharp interface is
`SharpInterface`) and co2lab's capillary fringe (`CapillaryFringeTable`:
`'P-scaled table'` and a corrected `'P-K-scaled table'`, on adaptive cubic
Hermite tables), with reconstruction, the inventory split and dissolution,
gated against a patched co2lab on `exampleVEJohansen.m`; and co2lab's
`'sharp_interface_integrated'` for rock varying with depth
(`SharpInterfaceIntegrated`), gated on `sharpInterfaceCrossSectionExample.m`;
and, not in co2lab, the fringe and the layers together (`LayeredFringe`: a
per-cell sum over each column's fine cells, built on integrals over any depth
interval for Stage 12), checked against both in their limits and against
quadrature. `examples/johansen_ve.jl` runs `exampleVEJohansen.m`.
Readable walkthroughs: `examples/sloping_aquifer.jl` (Stage 2),
`examples/residual_trapping.jl` (Stage 3, residual trapping) and
`examples/dissolution.jl` (Stage 4, rate-driven and instantaneous
dissolution), and `examples/atlas/` (Stage 5, the CO2 Storage Atlas with
plots: `show_co2_atlas.jl`, `models_from_atlas.jl`,
`show_top_surface_grids.jl`, `sleipner_model.jl`) and `examples/trapping/`
(Stage 8, co2lab-spillpoint's seven static examples and, with Stage 11's
interactive trap viewer, `showTrapsInteractively`) and
`examples/inventory/` (Stage 9, the trapping inventory of a run);
`examples/demos/` (Stage 11, interactive demos: the trap viewer,
co2lab's `exploreCapacity` and `exploreSimulation`);
`examples/layered_fringe.jl` (Stage 6, the capillary fringe in layered rock,
`LayeredFringe`, against the single-effect models and 3D) and
`examples/hybrid_grid.jl` (Stage 12, fine 3D cells around the injector, VE
elsewhere, against pure VE and 3D). The plan
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
| `select_secondary_variables!` and `@jutul_secondary` | saturations, `TotalMasses`, dissolved CO2; plume heights and relperm from the vertical model (`PlumeHeights`, `VERelPerm`; one method per model type, so the sharp interface reads no densities), `pcWG`, linear PVT, pore volume | `co2lab_fluid.jl`, `ve_fluid.jl`, `rock_compressibility.jl` |
| parameters (`default_parameter_values`) | `ColumnHeight`, `Rugosity` from the domain; without hysteresis `sGmax`, updated after each converged step (`update_after_step!`) | `ve_fluid.jl`, `co2lab_fluid.jl` |
| extra equations | `DissolvedCO2Equation` (rate-driven dissolution) and `HysteresisEquation` (`sGmax` with hysteresis), each with its own convergence criterion | `co2lab_fluid.jl` |
| `update_primary_variable!` | co2lab's switching unknown for instantaneous dissolution (`DissolvedOrFreeCO2`, on JutulDarcy's `BlackOilX`) | `instant_dissolution.jl` |
| cross term | `DissolutionWellSupplyCT`: a well's CO2 as dissolution supply in its cells | `co2lab_fluid.jl` |
| `apply_forces_to_equation!` | `VEBoundaryCondition` passed as JutulDarcy's `bc` force | `vebc.jl` |
| flux | none: JutulDarcy's stock two-point flux with the caprock `:gdz` is co2lab's VE flux (roadmap Stage 0) | — |
| wells | none: JutulDarcy's wells with co2lab's well index and a CO2 wellbore column (`co2_injector_control`) | `vewells.jl` |

The VE variables branch per cell on `is_ve`, and fluxes will branch per face
through a custom `kgrad` type: these are where Stage 12 (hybrid) adds fine 3D
cells (`../analysis/ve-meshing/isve-branch-proposal.md`).

Deliberate differences from co2lab are listed in the roadmap's standing
decisions ("Deliberate deviations from co2lab").

## Environment

`Jutul` and `JutulDarcy` are `dev`-ed from sibling checkouts, not from the
registry — the port tracks those working trees, and several roadmap decisions
cite them by `file:line`:

- `Jutul`: `../../Jutul.jl-fork`, the researcher's fork (github.com/Oddan/Jutul.jl)
  checked out on branch `co2lab-jutul-merge` (`aa88f180`) = the pinned upstream
  commit `36166ba6` (Jutul 0.4.31) plus the Jutul fixes this project needs: the
  cross-term row stride (#358) and signed tetrahedron volumes (#360). Both are
  merged upstream, but only in the Jutul 0.5 line, which no JutulDarcy release
  accepts yet (2026-10-06). `../../Jutul.jl` (upstream) is kept at `36166ba6`
  for reference.
- `JutulDarcy`: `../../JutulDarcy.jl`, the release 0.3.11 (`2c92d4d5`),
  unmodified.
- `GeoEnergyIO` (ECLIPSE input and corner-point grids, used through
  JutulDarcy): the registered release in `Manifest.toml` (1.1.34), not
  dev-ed. Up to 1.1.36 it tags some boundary faces with the wrong side (those
  facing an inactive cell of an active neighbouring column); JutulCO2Lab
  corrects them where it builds such grids (`johansen_ve_domain`). The fix is
  sintefmath/GeoEnergyIO.jl PR #66, from the researcher's fork
  (github.com/Oddan/GeoEnergyIO.jl, branch `boundary-face-direction`, cloned
  at `../../GeoEnergyIO.jl-fork`); the fork is only for that PR, nothing here
  uses it.

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

### Installing on another machine

1. Put the checkouts side by side (the manifests use these relative paths):

   ```
   <parent>/
     co2lab-jutul-merge-work/   this project
     Jutul.jl-fork/             Jutul, the fork's working branch
     JutulDarcy.jl/             JutulDarcy 0.3.11
     mrst-data/                 data sets for most examples (step 4)
     mrst/                      only for the gates against MRST (step 6)
     GeoEnergyIO.jl-fork/       not needed (only for GeoEnergyIO PR #66)
   ```

2. Clone them:

   ```
   cd <parent>
   git clone https://github.com/Oddan/co2lab-jutul-merge-work.git   # private
   git clone -b co2lab-jutul-merge https://github.com/Oddan/Jutul.jl.git Jutul.jl-fork
   git clone https://github.com/sintefmath/JutulDarcy.jl.git JutulDarcy.jl
   git -C JutulDarcy.jl checkout v0.3.11
   ```

3. With Julia 1.12 (1.12.6 here), instantiate the package's and the examples'
   environments, then run the tests (6–7 minutes on 8 threads):

   ```
   cd <parent>/co2lab-jutul-merge-work
   julia --project=JutulCO2Lab -e 'using Pkg; Pkg.instantiate()'
   julia --project=JutulCO2Lab/examples -e 'using Pkg; Pkg.instantiate()'
   julia --threads=8 --project=JutulCO2Lab -e 'using Pkg; Pkg.test()'
   ```

   `scripts/julia` pins Julia to the researcher's laptop's P-cores with
   `taskset`; elsewhere call `julia --threads=N` directly (or set
   `JULIA_CORES`). Tests needing data sets (step 4) or MRST reference dumps
   (step 6) are skipped until those are in place.

4. Data sets, in `<parent>/mrst-data/` (or a directory named by the
   environment variable `CO2LAB_ATLAS_DATA`), laid out as MRST's data
   directory:

   - `Johansen/`: `NPD5.grdecl`, `NPD5_Permeability.txt`,
     `NPD5_Porosity.txt` (the sector model, read directly by
     `johansen_ve_domain`) and `FULLFIELD_IMAXJMAX.GRDECL` with
     `FULLFIELD_FAULTS_IMAXJMAX.INC` (the full-field trap analysis). From
     MRST's Johansen data set: `mrst/core/examples/data/Johansen/` in an MRST
     checkout, or
     https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/johansen.zip.
   - `CO2Atlas/`, `CO2AtlasNorwegianSea/`, `CO2AtlasBarentsSea/`,
     `IGEMS surfaces/`, `IGEMSsample/`: the CO2 Storage Atlas and IGEMS
     surfaces, as MRST's dataset manager stores them (copy them from a
     machine that has them).
   - `Sleipner/`: the IEAGHG Sleipner benchmark (no longer downloadable; a
     copy from the researcher).

   A missing data set gives an error naming the directory (and, for
   Johansen, the download address). Terms of use of the data sets are
   unclear: copy them, do not publish them.

5. Examples: `julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/<file>.jl`
   (`JUTULCO2LAB_BACKEND=CairoMakie` writes figures to files). Without any
   data sets: `sloping_aquifer.jl`, `residual_trapping.jl`, `dissolution.jl`
   (co2lab's sloping aquifer, built by `sloping_aquifer`), `layered_fringe.jl`,
   `hybrid_grid.jl`, `atlas/show_top_surface_grids.jl`,
   `trapping/first_trapping_example.jl`. With `Johansen/`: `johansen_ve.jl`,
   `trapping/traps_johansen.jl`. The other `atlas/`, `trapping/`,
   `inventory/` and `demos/` examples need the Atlas, IGEMS or Sleipner data.

6. The gates against MRST (optional) need Octave (11.1 here; not MATLAB, whose
   Atlas grids differ), MRST on `gigaccs-working-branch`
   (`git clone -b gigaccs-working-branch https://github.com/SINTEF-AppliedCompSci/MRST.git mrst`,
   reference commit `50d892ad`) and `mrst-data/` (step 4). The reference dumps
   (`analysis/reference-data/*.mat`) are not in git: the gate scripts
   regenerate missing ones through Octave, the `scripts/mrst_dump_*.m`
   scripts the rest.

Once a JutulDarcy release accepts Jutul 0.5 and Jutul 0.5.2 (with #360) is
released, the fork can give way to registered versions.

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
  ve_fluid.jl        co2lab's VE fluid: vertical models (sharp interface), plume
                     heights, relperm, pcWG, linear PVT; column height and
                     rugosity parameters
  capillary_fringe.jl  the capillary-fringe vertical model: adaptive tables,
                     capillary pressure from saturation, relperm, fine-scale
                     profile, touched brine
  sharp_interface_integrated.jl  the integrated sharp interface: heights by
                     pore volume, relperm by permeability, over the fine cells
  layered_fringe.jl  the fringe in layered columns: per-cell sums over each
                     column's fine cells, integrals over any depth interval
  fine_hysteresis.jl  co2lab's endpoint-scaling hysteresis for fine 3D cells
                     (JutulDarcy relperm hysteresis and capillary pressure)
  fine_fluid.jl      the fluid at the fine scale: rock types' curves, the
                     hysteresis law; the VE fluid built from it; hybrid
                     grids' fine cells
  vehybrid.jl        hybrid grids' faces: categories, virtual-cell intervals,
                     transmissibilities
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
  ve_top_surface.jl  the caprock in map view (column polygons), contour lines

src/atlas/           the CO2 Storage Atlas (Stage 5)
  atlas.jl           data location, map reader, formation names, average rock,
                     coastline and wells
  projection.jl      UTM zone conversion (coastline to the Atlas's zone)
  sleipner.jl        the Sleipner benchmark (IEAGHG, layer 9) as a VE model
  johansen.jl        the Johansen NPD models (sector, full field), cut as co2lab
  atlas_grid.jl      top and thickness maps to a corner-point grid section, rock maps
  atlas_domain.jl    a formation as a VE domain: mesh, largest piece, columns,
                     rock; subsets of grid sections

src/plotting.jl      plotting functions (declared and documented; implemented in
                     the Makie extension)
ext/JutulCO2LabMakieExt/  Makie extension: VE map views, Atlas maps, coastline,
                     trap maps, in 2D or on the caprock in 3D (ve_axis3), and
                     the interactive trap viewer (interactive_traps, Stage 11:
                     GLMakie window or WGLMakie browser page); loaded with any
                     Makie backend

src/trapping/        structural traps (Stage 8, co2lab-spillpoint)
  trap_surface.jl    the caprock as trap analysis sees it (co2lab's Gt quantities),
                     stitched across faults
  cell_traps.jl      the cell-based method: traps, spill depths, rivers, regions
  node_traps.jl      the node-based method: spill field, regions, cycle merging
  node_to_cells.jl   node-based traps projected to columns (co2lab's n2cTraps)
  trap_analysis.jl   the front end, trap_analysis (co2lab's trapAnalysis)
  trap_quantities.jl trap volumes, downstream traps, capacity trees, injection
                     point, flattened caprock
  surfaces.jl        IRAP surfaces, IGEMS grids, grids under lattice surfaces
  capacity.jl        the storage capacity per column (co2lab-demo's exploreCapacity
                     estimate), the reachable structural capacity

src/ve/injection_case.jl  co2lab-demo's exploreSimulation scenario: wells, closed,
                     semi-open and open boundary stretches, a run handing over
                     each report step; the boundary in loops (Stage 11)

src/inventory/       where the CO2 is (Stage 9)
  trapping_inventory.jl  co2lab's trapping inventory per state and over a run

src/precompile.jl    precompiled workload: the trap analysis on a small synthetic
                     caprock (no data files)

examples/            readable walkthroughs, with their own environment
                     (examples/Project.toml: JutulCO2Lab, the Jutul fork,
                     GLMakie, CairoMakie), run from the project root with
                     scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/...
  atlas/             CO2 Storage Atlas scripts; JUTULCO2LAB_BACKEND=CairoMakie
                     writes figures to ./figures instead of opening windows
  trapping/          co2lab-spillpoint's examples (and showTrapsInteractively,
                     the interactive trap viewer)
  inventory/         a run's trapping inventory
  demos/             interactive demos (Stage 11), in a GLMakie window or a
                     browser (JUTULCO2LAB_BACKEND=WGLMakie): trap_viewer.jl (the
                     Atlas formations in the trap viewer), explore_capacity.jl
                     (co2lab's exploreCapacity), explore_simulation.jl
                     (co2lab's exploreSimulation)
test/runtests.jl     Oracle, Stage 1 (vs MRST dumps), fluid system, VE fluid, wells,
                     boundary conditions and reconstruction (vs MRST), Atlas, trapping,
                     inventory,
                     toolchain smoke test
test/atlas_checks.jl the Atlas comparisons with co2lab, shared with the Stage 5 gate
test/trap_checks.jl  the trapping comparisons with co2lab, shared with the Stage 8 gate
test/inventory_checks.jl  the inventory comparisons with co2lab, shared with the
                     Stage 9 item 2 check
../scripts/          MRST/Octave reference dumps, prototypes, gate scripts
```

All types and returned objects, in order of dependency: `STRUCTURES.md`.

## Open decisions

Tracked in `../analysis/STATUS.md`.

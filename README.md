# JutulCO2Lab

Vertical-equilibrium (VE) CO2 storage simulation on
[Jutul](https://github.com/sintefmath/Jutul.jl) /
[JutulDarcy](https://github.com/sintefmath/JutulDarcy.jl): a reimplementation
of MRST-co2lab, verified against it.

What it covers:

- **VE models on 3D grids**: columns (vertically connected stretches of a
  pillar, stackable), thickness-averaged rock, the caprock datum, boundary
  sides; co2lab's fluid system, VE relperm, capillary pressure, PVT (linear,
  or JutulDarcy's CO2–brine tables per column temperature), wells, pressure
  boundary conditions and reconstruction of 3D states.
- **Vertical models**: co2lab's sharp interface, its capillary fringe
  (`CapillaryFringeTable`: `'P-scaled table'` and a corrected
  `'P-K-scaled table'`), the integrated sharp interface for rock varying with
  depth (`SharpInterfaceIntegrated`), and, not in co2lab, the fringe in
  layered rock (`LayeredFringe`).
- **Residual trapping and dissolution**: co2lab's endpoint scaling and a
  capped fixed residual (`FixedResidual`); rate-driven and instantaneous
  dissolution.
- **Hybrid grids**: fine 3D cells among VE columns, with a transition flux.
- **The CO2 Storage Atlas, the Johansen models and the Sleipner benchmark**
  as VE models.
- **co2lab-spillpoint's trapping analysis** (both methods, per caprock
  sheet), the **trapping inventory** of a run, the **storage capacity**
  estimate, and interactive demos (trap viewer, co2lab-demo's
  `exploreCapacity` and `exploreSimulation`).

Every function's docstring states what it computes and returns; all types and
returned objects are listed in dependency order in `STRUCTURES.md`.

## Installation

JutulCO2Lab currently needs a Jutul with fixes not yet in a Jutul release
accepted by JutulDarcy, and a pinned JutulDarcy and GeoEnergyIO. Put the
checkouts side by side:

```
<parent>/
  JutulCO2Lab/            this repository
  Jutul.jl-fork/          Jutul, branch co2lab-jutul-merge
  JutulDarcy.jl/          JutulDarcy at ef1cfca5 (0.3.12 + two fixes)
  GeoEnergyIO.jl-fork/    GeoEnergyIO, branch co2lab-jutul-merge
  mrst-data/              data sets for most examples (below)
```

```
cd <parent>
git clone https://github.com/sintefmath/JutulCO2Lab.git
git clone -b co2lab-jutul-merge https://github.com/Oddan/Jutul.jl.git Jutul.jl-fork
git clone -b co2lab-jutul-merge https://github.com/Oddan/GeoEnergyIO.jl.git GeoEnergyIO.jl-fork
git clone https://github.com/sintefmath/JutulDarcy.jl.git JutulDarcy.jl
git -C JutulDarcy.jl checkout ef1cfca5
```

Then, with Julia 1.12, use the three checkouts in the package's and the
examples' environments, and run the tests (several minutes):

```
cd <parent>/JutulCO2Lab
julia --project=. -e 'using Pkg; Pkg.develop([PackageSpec(path = p) for p in
    ("../Jutul.jl-fork", "../JutulDarcy.jl", "../GeoEnergyIO.jl-fork")])'
julia --project=examples -e 'using Pkg; Pkg.develop([PackageSpec(path = p) for p in
    (".", "../Jutul.jl-fork", "../JutulDarcy.jl", "../GeoEnergyIO.jl-fork")])'
julia --threads=8 --project=. -e 'using Pkg; Pkg.test()'
```

`Manifest.toml` is not versioned: each checkout makes its own. The test suite
ends with a toolchain smoke test (a small stock immiscible case); if that
fails, the environment is the problem. Tests comparing with MRST reference
dumps or needing data sets are skipped when those are missing.

What the forks carry: Jutul 0.4.31 (the version JutulDarcy 0.3.12 accepts)
plus fixes merged upstream in the Jutul 0.5 line (#358, #360, #380) and
strided threading of secondary variables (#385); GeoEnergyIO 1.1.34 plus a
boundary-face fix (#66, in 1.2.0). Once a JutulDarcy release accepts Jutul 0.5,
registered versions replace them.

### Data sets

In `<parent>/mrst-data/` (or the directory named by `CO2LAB_ATLAS_DATA`),
laid out as MRST's data directory:

- `Johansen/`: `NPD5.grdecl`, `NPD5_Permeability.txt`, `NPD5_Porosity.txt`
  (the sector model) and `FULLFIELD_IMAXJMAX.GRDECL` with
  `FULLFIELD_FAULTS_IMAXJMAX.INC` (the full-field trap analysis); from MRST's
  Johansen data set, e.g.
  https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/johansen.zip.
- `CO2Atlas/`, `CO2AtlasNorwegianSea/`, `CO2AtlasBarentsSea/`,
  `IGEMS surfaces/`, `IGEMSsample/`: the CO2 Storage Atlas and IGEMS
  surfaces, as MRST's dataset manager stores them.
- `Sleipner/`: the IEAGHG Sleipner benchmark (no longer downloadable).

A missing data set gives an error naming the directory. Terms of use of the
data sets are unclear: copy them, do not publish them.

## Examples

```
julia --project=examples examples/<file>.jl
```

GLMakie opens a window per figure; `JUTULCO2LAB_BACKEND=CairoMakie` writes the
figures to `./figures` instead (`WGLMakie`: a browser page, for the demos).

- Without data sets: `sloping_aquifer.jl` (co2lab's sloping aquifer, step by
  step), `residual_trapping.jl`, `dissolution.jl`, `layered_fringe.jl`
  (`LayeredFringe` against the single-effect models and 3D), `hybrid_grid.jl`
  (fine cells around the injector, VE elsewhere), `atlas/show_top_surface_grids.jl`,
  `trapping/first_trapping_example.jl`, `trapping/sheets.jl`,
  `trapping/layered_traps.jl`.
- With `Johansen/`: `johansen_ve.jl` (co2lab's `exampleVEJohansen.m`),
  `trapping/traps_johansen.jl`.
- With the Atlas, IGEMS or Sleipner data: the other `atlas/` (`show_co2_atlas.jl`,
  `models_from_atlas.jl`, `sleipner_model.jl`), `trapping/` (co2lab-spillpoint's
  examples, `showTrapsInteractively`), `inventory/` and `demos/` scripts.

## Architecture

A VE model is a JutulDarcy reservoir model whose grid is a `CoarseMesh` of
VE columns over the 3D grid, and whose black-oil system carries co2lab's
formulation and VE physics as replaced variables. Nothing in Jutul or
JutulDarcy is patched.

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
| `select_secondary_variables!` and `@jutul_secondary` | saturations, `TotalMasses`, dissolved CO2; plume heights and relperm from the vertical model (`PlumeHeights`, `VERelPerm`), `pcWG`, PVT, pore volume | `co2lab_fluid.jl`, `ve_fluid.jl`, `rock_compressibility.jl` |
| parameters (`default_parameter_values`) | `ColumnHeight`, `Rugosity` from the domain; without hysteresis `sGmax`, updated after each converged step | `ve_fluid.jl`, `co2lab_fluid.jl` |
| extra equations | `DissolvedCO2Equation` (rate-driven dissolution) and `HysteresisEquation` (`sGmax` with hysteresis) | `co2lab_fluid.jl` |
| `update_primary_variable!` | co2lab's switching unknown for instantaneous dissolution (`DissolvedOrFreeCO2`); fine cells' saturation updates stopped at relperm kinks (`FineGasSaturation`) | `instant_dissolution.jl`, `co2lab_fluid.jl` |
| cross term | `DissolutionWellSupplyCT`: a well's CO2 as dissolution supply in its cells | `co2lab_fluid.jl` |
| `apply_forces_to_equation!` | `VEBoundaryCondition` passed as JutulDarcy's `bc` force | `vebc.jl` |
| flux | JutulDarcy's two-point flux with the caprock `:gdz`, which is co2lab's VE flux; on hybrid grids and stacked columns a transition flux (`VEFlux`) | `transition_flux.jl` |
| wells | JutulDarcy's wells with co2lab's well index and a CO2 wellbore column (`co2_injector_control`) | `vewells.jl` |

## Comparison with MRST-co2lab

The package was developed against MRST-co2lab: each capability is checked by a
numerical comparison with co2lab **on identical input**, per quantity with a
declared tolerance (`JutulCO2Lab.Oracle`):

```julia
using JutulCO2Lab.Oracle

ref = read_mrst_dump("stage1_johansen.mat")
rep = OracleReport("Johansen grid"; mrst = ref)
check!(rep, "column height", H, ref["height"]; tol = STAGE1_TOL[:height], perm = c2m)
display(rep)
passed(rep) || error("comparison failed")
```

The oracle never guesses orderings (pass an explicit `perm`, Julia index →
MRST index), and an empty report counts as a failure. Each report records the
commits of the Jutul, JutulDarcy and MRST checkouts (`source_provenance()`).

The MRST reference dumps and the Octave scripts that produce them are not
part of this repository; the tests that use them look for them in
`../analysis/reference-data/` and are skipped when they are missing. The
comparisons were made with MRST on branch `gigaccs-working-branch` (commit
`50d892ad`), under Octave.

## Layout

```
src/JutulCO2Lab.jl     module root
src/oracle.jl          Oracle: MRST dump reader and tolerance-checked comparison
src/provenance.jl      source commits and dirty state, recorded in every report

src/ve/                the VE simulator
  vegrid.jl            VE columns on the 3D mesh: partition, trimming, geometry
  vefaces.jl           face directions and column sides (topology)
  vedomain.jl          VE domain on a CoarseMesh: rock, caprock datum, sides
  vehybrid.jl          hybrid grids' faces: categories, intervals, transmissibilities
  transition_flux.jl   the flux across transition faces
  rock_compressibility.jl  co2lab's linear pore volume
  ve_fluid.jl          co2lab's VE fluid: vertical models, plume heights, relperm,
                       pcWG, PVT
  pvt_tables.jl        tabulated phase properties in pressure and temperature
  capillary_fringe.jl  the capillary fringe ('P-scaled table')
  sharp_interface_integrated.jl  the sharp interface in layered rock
  layered_fringe.jl    the fringe in layered rock (LayeredFringe)
  hysteresis_law.jl    residual-trapping laws (EndpointScaling, FixedResidual)
  fine_fluid.jl        the fluid at the fine scale; hybrid grids' fine cells
  fine_hysteresis.jl   hysteresis in fine 3D cells
  co2lab_fluid.jl      co2lab's fluid system, unknowns, equations; model, state
                       and simulator setup
  instant_dissolution.jl  instantaneous dissolution
  vewells.jl           VE wells
  vebc.jl              VE boundary conditions
  vereconstruct.jl     3D states from VE states
  ve_top_surface.jl    the caprock in map view, contour lines
  injection_case.jl    co2lab-demo's exploreSimulation scenario
  sloping_aquifer.jl   co2lab's sloping aquifer case

src/atlas/             the CO2 Storage Atlas, Johansen and Sleipner
src/trapping/          structural traps (co2lab-spillpoint), sheets, capacity
src/inventory/         the trapping inventory
src/plotting.jl        plotting functions (implemented in the Makie extension)
ext/JutulCO2LabMakieExt/  Makie extension: map views, Atlas maps, trap maps,
                       the interactive trap viewer
examples/              walkthroughs, with their own environment
test/                  the test suite (runtests.jl and the shared *_checks.jl)
```

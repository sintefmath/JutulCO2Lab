# Structures and objects in JutulCO2Lab

Every type, type alias, constant and returned composite object defined in the
package, in order of dependency: an entry depends only on entries above it
(built from them, holding them, or — for Jutul variables — computed from
them). Jutul/JutulDarcy supertypes are in brackets. Exported names are in
**bold**; the others are internal.

Keep this file current: update it in the same commit as any change that adds,
removes, renames or changes the meaning of a structure (see `../CLAUDE.md`).

The *File* column gives file names only: sections A–E are in `src/ve/`,
G in `src/atlas/`, F in `src/`.

*Last updated: 2026-09-26.*

## A. Grid and domain

| # | Name | File | What it is | Built from / depends on |
|---|---|---|---|---|
| 1 | **`VEColumns`** | `vegrid.jl` | Column partition of the 3D mesh: `cells` (per column, top to bottom), `pos`, `partition` (fine cell → column), `pillar` (logical `(i, j)`; not unique when stacked). A column is a vertically connected stretch of one pillar — in pure VE usually the whole pillar; pillars are cut where consecutive cells share no (unsealed) face | the 3D mesh (`ve_columns`, with `below_holes = :discard` or `:stack`) |
| 2 | *column geometry* (NamedTuple) | `vegrid.jl` | `(z_top, H, dz, z, dx, dy, z_cell)`: caprock depth, height, lateral extent per column; thickness, cumulative depth, centroid depth per fine cell | mesh, `VEColumns` (`ve_column_geometry`) |
| 3 | *averaged rock* (NamedTuple) | `vedomain.jl` | `(perm, poro, ntg)` per column, thickness-weighted; perm from x only | `VEColumns`, `dz` of 2 (`ve_average_rock`) |
| 4 | *VE domain* (NamedTuple) | `vedomain.jl` | `(domain, coarse, geometry, rock)`: the `DataDomain` on a `CoarseMesh` (with `:column_height`, `:rugosity`, `:is_ve` (all true until Stage 8), `:pillar_columns` (columns in the pillar) per column, caprock `:gdz`; `:transmissibilities` when stacked, sealing faces between stacked columns), that `CoarseMesh`, 2 and 3 | mesh, `VEColumns`, 2, 3, caprock rugosity (`ve_domain`) |
| 5 | **`VEColumnSides`** | `vefaces.jl` | Each column's lateral fine faces grouped by direction and neighbouring column: one entry per (column, direction, neighbour) — one coarse face seen from one side — or an outer side (neighbour 0) when the side has no neighbour; `exposed_faces`: boundary fine faces on sides that have a neighbour (no-flow) | mesh (direction tags), `VEColumns` (`ve_column_sides`) |
| 6 | **`VEBoundarySides`** | `vedomain.jl` | Outer sides (no neighbour): direction, caprock edge depth and midpoint, area, centroid, normal, half transmissibility, fine faces, `across_fault` | 5, `VEColumns`, 4's domain (`ve_boundary_sides`, which also returns the unassigned faces) |
| 7 | **`VEBoundaryCondition`** [force, passed as JutulDarcy's `bc`] | `vebc.jl` | co2lab's pressure BC on one open side: column, boundary pressure, transmissibility, edge-minus-caprock depth `dz`, inflow saturations, `carry_dissolved` (with dissolution: brine carries its dissolved CO2 across, or co2lab's booking as brine) | 6 (selected with `ve_select_sides`), 2 (`ve_boundary_conditions`) |

Constants: `_DIRECTIONS` (`:left, :right, :lower, :upper, :top, :bottom` —
Jutul's direction names in MRST's face-tag order), `_TOP`, `_BOTTOM`
(`vefaces.jl`).

## B. Fluid settings

| # | Name | File | What it is | Built from / depends on |
|---|---|---|---|---|
| 8 | **`CO2LabVEFluid`** | `ve_fluid.jl` | co2lab's VE fluid settings (makeVEFluid's fluid options): residual saturations, endpoints, reference densities and linear PVT, viscosities, dissolution (`dissolution`, `dis_rate`, `dis_max`; `dis_rate` 0 or `Inf` is instantaneous). Rugosity is in 4, rock compressibility in the model (20) | nothing (`co2lab_ve_fluid`); supplies the constants of the variables in C |

## C. Jutul variables and equations, in the order they are computed

| # | Name | File | Role | Computed from |
|---|---|---|---|---|
| 9 | `ColumnHeight` [ScalarVariable] | `ve_fluid.jl` | parameter `H` | `domain[:column_height]` (4) |
| 9b | `Rugosity` [ScalarVariable] | `ve_fluid.jl` | parameter `dh`, caprock rugosity per column | `domain[:rugosity]` (4); zero if absent |
| 10 | **`MaxGasSaturation`** [ScalarVariable] | `co2lab_fluid.jl` | `sGmax` in `[0, 1 − res_water]`: without hysteresis a parameter, the running maximum of converged `sG`; with hysteresis (`res_gas > 0`) an unknown, governed by 22c | `sG`; constant `res_water` from 8 |
| 11 | `GasSaturation` [ScalarVariable] | `co2lab_fluid.jl` | unknown `sG` | — |
| 12 | `DissolvedGasRatio` [ScalarVariable] | `co2lab_fluid.jl` | unknown `rs` (rate-driven dissolution) | — |
| 12b | **`DissolvedOrFreeCO2`** [ScalarVariable, values JutulDarcy's `BlackOilX`] | `instant_dissolution.jl` | co2lab's switching unknown `:CO2X` for instantaneous dissolution: `rs` where brine is undersaturated, `sG` where saturated; co2lab's update rule | — |
| 12c | `InstantGasSaturation`, `InstantDissolvedGasRatio` [ScalarVariable] | `instant_dissolution.jl` | secondary `sG`, `Rs` from 12b (instantaneous dissolution) | 12b |
| 13 | `PlumeHeights` [VectorVariables] | `ve_fluid.jl` | secondary: row 1 `h`, row 2 `h_max`; field `is_ve` (`nothing` = all VE; 0 in fine cells) | 11, 10, 9 |
| 14 | `VESharpInterfaceRelPerm` [AbstractRelativePermeabilities] | `ve_fluid.jl` | relative permeabilities (co2lab's `sharp_interface_simple`); fields `is_ve` and `fine` (3D relperm for fine cells, Stage 8) | 13, 9, 9b |
| 15 | **`LinearShrinkage`** [PhaseVariables] | `co2lab_fluid.jl` | shrinkage factors `b` | pressure |
| 16 | `BlackOilDensities` [PhaseVariables] | `co2lab_fluid.jl` | phase mass densities | 15 (and 12) |
| 17 | `VECapillaryPressure` [AbstractCapillaryPressure] | `ve_fluid.jl` | `pcWG = Δρ·g·h`; field `is_ve` (fine cells: Stage 8) | 13, 16 |
| 18 | `LinearViscosities` [PhaseVariables] | `ve_fluid.jl` | viscosities (used by `setup_co2lab_model`) | pressure |
| 19 | **`ConstantViscosities`** [PhaseVariables] | `co2lab_fluid.jl` | placeholder default, replaced by 18 in `setup_co2lab_model` | — |
| 20 | `LinearPoreVolume` [ScalarVariable] | `rock_compressibility.jl` | pore volume `1 + c·(p − p_ref)` times the static pore volume — co2lab's rock compressibility, set on the model by `set_co2lab_rock_compressibility!` (or `setup_co2lab_model`'s `rock_compressibility`); exchangeable for JutulDarcy's `LinearlyCompressiblePoreVolume` | pressure, parameter `StaticFluidVolume` |
| 21 | `DissolvedCO2Mass` [ScalarVariable] | `co2lab_fluid.jl` | dissolved CO2 mass per cell | 12, 15, saturations, 20 |
| 22 | **`DissolvedCO2Equation`** [JutulEquation] | `co2lab_fluid.jl` | dissolved-CO2 conservation with co2lab's transfer (`dis_rate`, `dis_max`; demand and supply caps, touched-brine minimum) — rate-driven dissolution | 21, fluxes, 13 |
| 22b | `InstantDissolution` (marker) | `instant_dissolution.jl` | tells 22c the depletion is the dissolved-CO2 balance's residual (instantaneous dissolution) | — |
| 22c | **`HysteresisEquation`** [JutulEquation] | `co2lab_fluid.jl` | co2lab's `sGmax − max(sGmax0 − depletion, sG) = 0` per cell (hysteresis on); `dissolution`: `nothing`, 22 (depletion = its transfer) or 22b | 10, 11, 22 or 22b |
| 22d | `DissolutionWellSupplyCT` [AdditiveCrossTerm] | `co2lab_fluid.jl` | well → reservoir cross term on 22 (and 22c): the well's CO2 inflow counted as dissolution supply in its perforated cells | 22, 22c, the well's state |

## D. System and model types

| # | Name | File | What it is | Depends on |
|---|---|---|---|---|
| 23 | `ConstantRs` | `co2lab_fluid.jl` | co2lab's constant saturated `rsSat` | — |
| 24 | `CO2LabSystem` (alias) | `co2lab_fluid.jl` | JutulDarcy's `StandardBlackOilSystem` with formulation tag `:co2lab`; the variables of C attach to it | 23 when dissolution is on |
| 25 | `CO2LabNoDissolutionSystem` (alias) | `co2lab_fluid.jl` | 24 without dissolved gas | 24 |
| 26 | `CO2LabModel`, `CO2LabNoDissolutionModel` (aliases) | `co2lab_fluid.jl` | `SimulationModel`s of 24 / 25, for method dispatch | 24, 25 |
| 27 | *co2lab model* (`MultiModel`) | `co2lab_fluid.jl` | `:Reservoir` (26, with C), one submodel per well, `:Facility` | 4, 8, rock compressibility, well domains (`setup_co2lab_model`) |

## E. Wells and results

| # | Name | File | What it is | Built from / depends on |
|---|---|---|---|---|
| 28 | *VE well* (JutulDarcy well `DataDomain`) | `vewells.jl` | simple well in one column: co2lab's well index, perforation at the caprock, node at the chosen datum | 4 (`setup_ve_well`, `ve_well_index`) |
| 29 | *CO2 injector control* (JutulDarcy `InjectorControl`) | `vewells.jl` | CO2 injection with the reference density and a CO2 wellbore column; also fills the wellbore at the start (`setup_co2lab_state(...; forces)`) | 8 (`co2_injector_control`) |
| 30 | *reconstructed 3D state* (NamedTuple) | `vereconstruct.jl` | `(pressure, sG, sGmax, rs)` per 3D cell; `rs` `nothing` without dissolution | a reservoir state of 27, 4, `VEColumns`, 8 (`ve_reconstruct`) |
| 31 | *boundary mass fluxes* (Matrix) | `vebc.jl` | 3 × number of conditions, kg/s out: brine, CO2 (free and dissolved), the dissolved part | 7, a reservoir state (`ve_boundary_mass_fluxes`) |

Constants: `_WELL_CONSTANT_RATIO`, `_WELL_CONSTANT_VALUE` (MRST's
`wellConstant` table, `vewells.jl`).

## G. CO2 Storage Atlas (roadmap Stage 5), independent of A–E

| # | Name | File | What it is | Built from / depends on |
|---|---|---|---|---|
| 37 | *ESRI ASCII grid* (NamedTuple) | `atlas.jl` | `(meta, data)`: header `(ncols, nrows, xllcorner, yllcorner, cellsize, nodata)`; `data` ncols × nrows, x then y from the south, NaN for missing/zero | an Atlas file (`read_aai_grid`; `co2lab_rows` reproduces co2lab's one-row shift) |
| 38 | *Atlas map file* (NamedTuple) | `atlas.jl` | `(name, variant, path)`: co2lab's naming (`Huginfm_east_thickness` → `Huginfmeast`, `:thickness` or `:top`) | the data directory (`atlas_map_files`) |
| 39 | *Atlas map* (NamedTuple) | `atlas.jl` | `(name, variant, meta, data)`: 37 with tops as depths, coarsened (and refined, `atlas_grid.jl`); `meta` adds `dims` and `x0`, `y0`, the position of `data[1, 1]` (its cell centre, or co2lab's corner with `co2lab_compatible`); `cellsize` is the spacing of the values | 38, 37 (`read_atlas_maps`) |
| 40 | *Atlas average rock* (NamedTuple) | `atlas.jl` | `(perm, poro, ntg)`: the Atlas reports' averages, perm in m², NaN where none, `ntg` `nothing` unless given | formation name (`atlas_average_rock`) |
| 41 | *Atlas grid section* (NamedTuple) | `atlas_grid.jl` | `(grid, petro)`: `grid` a GeoEnergyIO grid section `Dict{String, Any}` (`cartDims`, `COORD`, `ZCORN` with 0 at undefined pillars, `ACTNUM` `Array{Bool, 3}`, and `PORO`, `PERMX/Y/Z`, `NTG` as `Array{Float64, 3}` where the Atlas has rock maps); `petro` 40 with map means where there are maps | 39, 40 (`atlas_grid_section`) |

Constants: `_ATLAS_DATASETS` (the three data-set directories),
`_ATLAS_AVERAGE_ROCK` (co2lab's `getAvgRock` table), `_JURASSIC_MID`,
`_JURASSIC_TOP` (formations taking the `Jurassicmid` and `Jurassic` top maps,
`atlas_grid.jl`).

## F. MRST comparison harness (`JutulCO2Lab.Oracle`), independent of A–E

| # | Name | File | What it is | Depends on |
|---|---|---|---|---|
| 32 | **`MrstDump`** | `oracle.jl` | a loaded MRST export (`read_mrst_dump`) | — |
| 33 | **`SourceInfo`** | `provenance.jl` | commit and dirty state of one source tree | — |
| 34 | **`Tolerance`** | `oracle.jl` | relative and absolute tolerance; pass if `|a − b| ≤ atol + rtol·|b|` | — |
| 35 | `CheckResult` | `oracle.jl` | outcome of one comparison | 34 |
| 36 | **`OracleReport`** | `oracle.jl` | titled list of results plus provenance | 35, 33; records the MRST commit from 32 |

Constants: **`STAGE1_TOL`** (Stage 1 gate tolerances, 34), `TRACKED_PACKAGES`
(Jutul and JutulDarcy, for provenance).

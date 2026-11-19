"""
    JutulCO2Lab

MRST-co2lab's tools for CO2 storage, reimplemented on Jutul/JutulDarcy.

# Parts

- **VE simulation**: vertical-equilibrium (VE) grids and domains from 3D
  meshes ([`ve_columns`](@ref), [`ve_domain`](@ref)), also hybrid VE/3D
  grids; co2lab's fluid, relperm, capillary fringe, hysteresis and
  dissolution ([`co2lab_ve_fluid`](@ref), [`setup_co2lab_model`](@ref));
  wells, boundary conditions, and reconstruction of 3D states.
- **CO2 Storage Atlas**: the Norwegian Atlas formations, the Sleipner and
  Johansen models as meshes and VE domains ([`atlas_ve_domain`](@ref)).
- **Trapping analysis**: structural traps and spill paths of a caprock
  ([`trap_analysis`](@ref)) and storage capacity estimates.
- **Trapping inventory**: where the CO2 of a simulation is, by trapping
  mechanism ([`trapping_inventory`](@ref)).
- **Plotting**: maps, traps and inventories; implemented in the Makie
  extension, loaded with any Makie backend.
- **MRST comparison**: the `JutulCO2Lab.Oracle` harness, which checks
  results against MRST-co2lab on identical input with explicit tolerances.
"""
module JutulCO2Lab

using Jutul
using JutulDarcy
using LinearAlgebra
using SparseArrays
using Statistics
using Printf
using StaticArrays
using Graphs
using MAT
using PrecompileTools

# ------------------------------ MRST comparison -------------------------------

include("oracle.jl")
using .Oracle
export Oracle

# ---------------------------- VE grids and domains ----------------------------

include("ve/vegrid.jl")
export VEColumns, ve_mesh_from_mrst, ve_columns, remove_cells, ve_column_geometry
export ve_except_pillars
export number_of_columns, column_cells

include("ve/vefaces.jl")
export VEColumnSides, ve_column_sides

include("ve/vehybrid.jl")
include("ve/vedomain.jl")
export net_rock, ve_average_rock, ve_domain, VEBoundarySides, ve_boundary_sides
export ve_hydrostatic_pressure, ve_hydrostatic_profile, geothermal_gradient

# ----------------------------- VE fluid and model -----------------------------

include("ve/rock_compressibility.jl")
export set_co2lab_rock_compressibility!

include("ve/pvt_tables.jl")
export PhaseTable, jutuldarcy_property_tables

include("ve/ve_fluid.jl")
include("ve/capillary_fringe.jl")
include("ve/sharp_interface_integrated.jl")
include("ve/hysteresis_law.jl")
include("ve/layered_fringe.jl")
include("ve/fine_hysteresis.jl")
include("ve/fine_fluid.jl")
include("ve/transition_flux.jl")
export co2lab_ve_fluid, CO2LabVEFluid, instantaneous_dissolution, needs_temperature
export VerticalModel, SharpInterface, SharpInterfaceIntegrated, CapillaryFringeTable
export LayeredFringe, co2lab_capillary_curve
export EndpointScalingHysteresis, FreeSaturationCapillaryPressure
export HysteresisLaw, EndpointScaling, FixedResidual, RockCurves, FineScaleFluid

include("ve/co2lab_fluid.jl")
export co2lab_system, has_dissolution
export setup_co2lab_model, setup_co2lab_state, setup_co2lab_simulator
export LinearShrinkage, ConstantViscosities, MaxGasSaturation, DissolvedCO2Equation
export HysteresisEquation

include("ve/instant_dissolution.jl")
export DissolvedOrFreeCO2

# -------------------- Wells, boundaries, 3D states, cases ---------------------

include("ve/vewells.jl")
export ve_well_index, setup_ve_well, co2_injector_control

include("ve/vebc.jl")
export VEBoundaryCondition, ve_select_sides, ve_boundary_conditions
export ve_boundary_mass_fluxes

include("ve/vereconstruct.jl")
export ve_reconstruct

include("ve/injection_case.jl")
include("ve/sloping_aquifer.jl")
export setup_injection_case, simulate_injection_case, ve_boundary_loops
export sloping_aquifer_grid, sloping_aquifer, rampup_timesteps

include("ve/ve_top_surface.jl")
include("ve/ve_section.jl")
export ve_top_surface, ve_contour_lines, ve_section

# ----------------------------- CO2 Storage Atlas ------------------------------

include("atlas/projection.jl")
include("atlas/atlas.jl")
export atlas_data_directory, atlas_formation_names, read_aai_grid, atlas_map_files
export read_atlas_maps, atlas_average_rock, atlas_coastline, atlas_wells

include("atlas/atlas_grid.jl")
export atlas_grid_section, atlas_grid_formations, atlas_map_values

include("atlas/atlas_domain.jl")
export atlas_ve_domain, atlas_mesh, grid_section_subset

include("atlas/sleipner.jl")
export sleipner_data_directory, sleipner_grid_section, sleipner_ve_domain

include("atlas/johansen.jl")
export johansen_data_directory, johansen_mesh, johansen_well_point
export johansen_grid_section, johansen_ve_domain

# ----------------------------- Trapping analysis ------------------------------

include("trapping/trap_surface.jl")
include("trapping/sheets.jl")
include("trapping/cell_traps.jl")
include("trapping/node_traps.jl")
include("trapping/node_to_cells.jl")
include("trapping/trap_analysis.jl")
include("trapping/trap_quantities.jl")
include("trapping/surfaces.jl")
include("trapping/capacity.jl")
export trap_surface, trap_analysis, trap_volumes, downstream_traps
export ve_sheets, sheet_spill_targets
export maximize_trapping, optimal_injection_point, flatten_traps
export read_irap_surface, igems_surface_path, igems_mesh, surface_mesh
export trapping_capacity, reachable_capacity

# ----------------------------- Trapping inventory -----------------------------

include("inventory/trapping_inventory.jl")
export TRAPPING_CATEGORIES, trapping_inventory, trapping_history

# ------------------------- Plotting (Makie extension) -------------------------

include("plotting.jl")
export plot_ve_map, plot_ve_map!, plot_atlas_map, plot_atlas_map!
export plot_atlas_context!, plot_ve_outline!, plot_ve_surface!, show_figure
export plot_traps, plot_traps!, ve_axis3
export plot_trapping_history, plot_trapping_history!
export plot_ve_section, plot_ve_section!
export interactive_traps, interactive_traps!, show_column_values!

# ------------------------------- Precompilation -------------------------------

include("precompile.jl")

end # module

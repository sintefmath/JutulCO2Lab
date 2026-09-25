"""
    JutulCO2Lab

Reimplementation of MRST-co2lab's vertical-equilibrium (VE) CO2 storage
simulator on top of Jutul/JutulDarcy.

Exploratory. See `../analysis/roadmap.org` for the staged plan and
`../analysis/ve-meshing/hybrid-ve-feasibility.md` for the evidence base.

Contents: the oracle harness (`JutulCO2Lab.Oracle`) — every acceptance gate
in the roadmap is a numerical comparison against MRST-co2lab on identical
input; the VE grid and domain (Stage 1); co2lab's fluid system (black-oil
with a `:co2lab` formulation, Stage 2 item 0); co2lab's VE relperm, `pcWG`
and PVT (items 1, 3); VE wells (item 4); VE boundary conditions (item 5);
reconstruction of 3D states (item 6).
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

include("oracle.jl")
using .Oracle
export Oracle

include("vegrid.jl")
export VEColumns, ve_mesh_from_mrst, ve_columns, ve_trim, ve_column_geometry
export number_of_columns, column_cells

include("vefaces.jl")
export VEColumnSides, ve_column_sides

include("vedomain.jl")
export ve_average_rock, ve_domain, VEBoundarySides, ve_boundary_sides
export ve_hydrostatic_pressure

include("rock_compressibility.jl")
export set_co2lab_rock_compressibility!

include("ve_fluid.jl")
export co2lab_ve_fluid, CO2LabVEFluid

include("co2lab_fluid.jl")
export co2lab_system, has_dissolution
export setup_co2lab_model, setup_co2lab_state, setup_co2lab_simulator
export LinearShrinkage, ConstantViscosities, MaxGasSaturation, DissolvedCO2Equation
export HysteresisEquation

include("vewells.jl")
export ve_well_index, setup_ve_well, co2_injector_control

include("vebc.jl")
export VEBoundaryCondition, ve_select_sides, ve_boundary_conditions
export ve_boundary_mass_fluxes

include("vereconstruct.jl")
export ve_reconstruct

end # module

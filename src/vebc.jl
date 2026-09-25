# VE boundary conditions (roadmap Stage 2 item 5).
#
# co2lab's pressure boundary condition on the lateral sides of the VE columns,
# as MRST computes it (getBoundaryConditionFluxesAD.m, called by
# ReservoirModel.addBoundaryConditionsAndSources for CO2VEBlackOilTypeModel).
# Per phase, with the column's caprock pressure p (brine) and p + pcWG (CO2):
#
#   dP = p_bc − p_phase − ρ_phase·g·(z_edge − z_cell)
#
# - one boundary pressure p_bc for both phases (no capillary pressure outside);
# - the gravity term runs from the column's caprock depth to the side's
#   caprock-edge depth (MRST's Gt.cells.z and Gt.faces.z);
# - densities are the column's (MRST averages inside and outside, and the
#   outside density defaults to the inside one);
# - outflow (dP ≤ 0) with the phase's own mobility, inflow (dP > 0) with the
#   column's total mobility times the inflow saturation (brine only by
#   default: open_boundary_conditions.m's `'sat', [1, 0]`).
#
# JutulDarcy's FlowBoundaryCondition has none of this (no pcWG, no gravity
# term, one flow direction for all phases), hence a separate force. It is
# passed as JutulDarcy's `bc` force: setup_reservoir_forces(model; bc = bcs).

"""
    VEBoundaryCondition(cell, pressure, trans, dz, sat)

co2lab's pressure boundary condition on one lateral side of VE column `cell`:
boundary pressure `pressure` (Pa, both phases), transmissibility `trans` (m³,
permeability times geometry; the column's half transmissibility to the side),
`dz` = caprock-edge depth − column caprock depth (m, depth positive down),
`sat` = saturations (brine, CO2) of inflowing fluid. Build with
[`ve_boundary_conditions`](@ref).
"""
struct VEBoundaryCondition{T <: Real}
    cell::Int
    pressure::T
    trans::T
    dz::T
    sat::NTuple{2, T}
end

"""
    ve_select_sides(sides, columns; directions = (:left, :right, :lower, :upper),
                    i = nothing, j = nothing, include_across_fault = false)
        -> Vector{Int}

Indices into `sides` ([`ve_boundary_sides`](@ref)) of the sides to open, the
counterpart of MRST's `searchForBoundaryFaces`: sides in one of `directions`
(Jutul's direction names), whose column's logical `i` (and `j`) lies in `i`
(`j`) — a range, a collection or `nothing` for any. Sides flagged
`across_fault` (a fully offset fault, sealed) are left out unless
`include_across_fault = true`.

Returns a new `Vector{Int}`, increasing. `sides` and `columns` are only read.
"""
function ve_select_sides(sides::VEBoundarySides, columns::VEColumns;
                         directions = (:left, :right, :lower, :upper),
                         i = nothing, j = nothing, include_across_fault::Bool = false)
    dirs = map(d -> _direction_index(d), collect(directions))
    inrange(x, r) = isnothing(r) || x in r
    keep(s) = sides.direction[s] in dirs &&
              inrange(columns.pillar[sides.column[s]][1], i) &&
              inrange(columns.pillar[sides.column[s]][2], j) &&
              (include_across_fault || !sides.across_fault[s])
    return filter(keep, collect(1:length(sides)))
end

# 1..4 for a lateral direction name
function _direction_index(d::Symbol)
    k = findfirst(==(d), _DIRECTIONS)
    (isnothing(k) || k > 4) &&
        throw(ArgumentError("lateral direction must be one of " *
                            "$(_DIRECTIONS[1:4]), got :$d"))
    return k
end

"""
    ve_boundary_conditions(ve, sides, ix; pressure, sat = (1.0, 0.0),
                           trans = nothing) -> Vector{VEBoundaryCondition}

co2lab's pressure boundary conditions on the sides `sides[ix]` (from
[`ve_boundary_sides`](@ref) and [`ve_select_sides`](@ref)); `ve` is the
NamedTuple from [`ve_domain`](@ref).

- `pressure`: a function of the side's caprock-edge depth (e.g. hydrostatic,
  `z -> ρ_w·g·z`, as co2lab's `open_boundary_conditions`) or one value per
  selected side (Pa).
- `sat`: saturations (brine, CO2) of inflowing fluid; co2lab's default is
  brine.
- `trans`: one transmissibility per selected side; defaults to the column's
  half transmissibility to the side, `sides.half_trans`. Give MRST's
  `T_all(bc.face)` for runs on identical input.

Returns a new vector, one condition per selected side, in the order of
`ix`; pass it as `setup_reservoir_forces(model; bc = ...)`. Arguments are
only read.
"""
function ve_boundary_conditions(ve, sides::VEBoundarySides, ix; pressure,
                                sat = (1.0, 0.0), trans = nothing)
    n = length(ix)
    z_edge = sides.edge_depth[ix]
    p_bc = pressure isa Function ? pressure.(z_edge) : collect(Float64, pressure)
    T = isnothing(trans) ? sides.half_trans[ix] : collect(Float64, trans)
    length(p_bc) == n && length(T) == n ||
        throw(ArgumentError("need one pressure and one transmissibility per side"))
    s = (Float64(sat[1]), Float64(sat[2]))
    sum(s) ≈ 1 || throw(ArgumentError("inflow saturations must sum to 1, got $sat"))
    cols = sides.column[ix]
    dz = z_edge .- ve.geometry.z_top[cols]
    return [VEBoundaryCondition(cols[k], p_bc[k], T[k], dz[k], s) for k in 1:n]
end

"""
    _ve_bc_volume_flux(bc, p_ph, mob_ph, rho_ph, mob_t, sat_ph, g) -> flux

Reservoir-volume flux (m³/s) of one phase *into* the column through `bc`,
MRST's `qRes`: `dP = p_bc − p_ph − rho_ph·g·dz`; for `dP > 0` (inflow)
`mob_t·T·dP·sat_ph`, otherwise `mob_ph·T·dP`.
"""
function _ve_bc_volume_flux(bc, p_ph, mob_ph, rho_ph, mob_t, sat_ph, g)
    dP = bc.pressure - p_ph - rho_ph * g * bc.dz
    return dP > 0 ? mob_t * bc.trans * dP * sat_ph : mob_ph * bc.trans * dP
end

"""
    _ve_bc_mass_outflow(bc, state, g) -> (brine, co2)

Mass flux (kg/s) of brine and CO2 *out of* the column through `bc`, from the
column's pressure, capillary pressure, relative permeabilities, viscosities
and densities in `state` (read with `state[:Name]`, so both a simulator state
and an output state work).
"""
function _ve_bc_mass_outflow(bc, state, g)
    c = bc.cell
    kr, mu, rho = state[:RelativePermeabilities], state[:PhaseViscosities],
                  state[:PhaseMassDensities]
    p = state[:Pressure][c]
    # CO2 pressure = brine pressure + pcWG (JutulDarcy's non-reference phase)
    pp = (p, p + state[:CapillaryPressure][1, c])
    mob = (kr[1, c] / mu[1, c], kr[2, c] / mu[2, c])
    mob_t = mob[1] + mob[2]
    q(ph) = -rho[ph, c] *
            _ve_bc_volume_flux(bc, pp[ph], mob[ph], rho[ph, c], mob_t, bc.sat[ph], g)
    return (q(1), q(2))
end

"""
    ve_boundary_mass_fluxes(bcs, state) -> Matrix{Float64}

Mass fluxes out of the domain (kg/s, positive out) through each boundary
condition in `bcs`, evaluated on `state`, which must hold `:Pressure`,
`:CapillaryPressure`, `:RelativePermeabilities`, `:PhaseViscosities` and
`:PhaseMassDensities` (add them to the reservoir model's
`output_variables` to get them in simulation output). Dissolution off only.

Returns a new 2 × `length(bcs)` matrix: row 1 brine, row 2 CO2. Arguments
are only read.
"""
function ve_boundary_mass_fluxes(bcs::AbstractVector{<:VEBoundaryCondition}, state)
    q = zeros(2, length(bcs))
    for (k, bc) in enumerate(bcs)
        q_k = _ve_bc_mass_outflow(bc, state, Jutul.gravity_constant)
        q[1, k], q[2, k] = Jutul.value(q_k[1]), Jutul.value(q_k[2])
    end
    return q
end

# Modifies: `acc` (the mass-conservation accumulation, per component and cell).
# Components are the phases (dissolution off): brine, CO2.
function Jutul.apply_forces_to_equation!(
        acc, storage, model::CO2LabModel, eq::JutulDarcy.ConservationLaw{:TotalMasses},
        eq_s, force::AbstractVector{<:VEBoundaryCondition}, time)
    has_dissolution(model.system) &&
        error("VE boundary conditions with dissolution are roadmap Stage 5 " *
              "(co2lab passes no dissolved components to its boundary fluxes)")
    state = storage.state
    for bc in force
        q = _ve_bc_mass_outflow(bc, state, Jutul.gravity_constant)
        acc[1, bc.cell] += q[1]
        acc[2, bc.cell] += q[2]
    end
    return acc
end

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
#
# With dissolution, co2lab passes no dissolved components to the boundary
# fluxes (CO2VEBlackOilTypeModel.addBoundaryConditionsAndSources, dissolved =
# {}): the brine phase's mass — dissolved CO2 included — is booked to the
# brine component, and the dissolved-CO2 law sees no boundary flux — the
# brine balance loses the dissolved CO2's mass (analysis/co2lab-issues.org).
# Here brine leaves with its dissolved CO2 and enters without any
# (`dissolved = :carry`, the default); `dissolved = :co2lab` reproduces
# co2lab's booking for comparisons.

# ---------------- Boundary conditions and the sides they open -----------------

"""
    VEBoundaryCondition(cell, pressure, trans, dz, sat, carry_dissolved)

co2lab's pressure boundary condition on one lateral side of VE column `cell`:
boundary pressure `pressure` (Pa, both phases), transmissibility `trans` (m³,
permeability times geometry; the column's half transmissibility to the side),
`dz` = caprock-edge depth − column caprock depth (m, depth positive down),
`sat` = saturations (brine, CO2) of inflowing fluid, `carry_dissolved`: with
dissolution, whether brine carries its dissolved CO2 across (see the top of
`vebc.jl`). Build with [`ve_boundary_conditions`](@ref).
"""
struct VEBoundaryCondition{T <: Real}
    cell::Int
    pressure::T
    trans::T
    dz::T
    sat::NTuple{2, T}
    carry_dissolved::Bool
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
                           trans = nothing, dissolved = :carry)
        -> Vector{VEBoundaryCondition}

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
- `dissolved` (with dissolution only): `:carry` — brine leaves with its
  dissolved CO2 and enters without any; `:co2lab` — the brine phase's mass,
  dissolved CO2 included, is booked to brine, as co2lab does (for
  comparisons with co2lab).

Returns a new vector, one condition per selected side, in the order of
`ix`; pass it as `setup_reservoir_forces(model; bc = ...)`. Arguments are
only read.
"""
function ve_boundary_conditions(ve, sides::VEBoundarySides, ix; pressure,
                                sat = (1.0, 0.0), trans = nothing, dissolved = :carry)
    dissolved in (:co2lab, :carry) ||
        throw(ArgumentError("dissolved must be :co2lab or :carry, got $dissolved"))
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
    carry = dissolved == :carry
    return [VEBoundaryCondition(cols[k], p_bc[k], T[k], dz[k], s, carry) for k in 1:n]
end

# ----------------------------------- Fluxes -----------------------------------

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
    _ve_bc_mass_outflow(bc, state, g) -> (brine, co2, dissolved)

Mass flux (kg/s) *out of* the column through `bc`, from the column's
pressure, capillary pressure, relative permeabilities, viscosities and
densities in `state` (read with `state[:Name]`, so both a simulator state and
an output state work), per component row of the mass balance: brine and CO2,
and the part of the CO2 that is dissolved (zero without dissolution, and with
co2lab's booking, where the brine phase's whole mass counts as brine).
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

    vol(ph) = -_ve_bc_volume_flux(bc, pp[ph], mob[ph], rho[ph, c], mob_t, bc.sat[ph], g)
    brine, co2 = rho[1, c] * vol(1), rho[2, c] * vol(2)
    bc.carry_dissolved && haskey(state, :Rs) || return (brine, co2, zero(brine))

    # brine leaves with its dissolved CO2 and enters without any; with
    # ρ_L = b_L·(ρ_LS + rs·ρ_VS) and ρ_V = b_V·ρ_VS, the dissolved part of the
    # brine phase's mass is b_L·rs·ρ_VS per volume
    b, rs = state[:ShrinkageFactors], state[:Rs][c]
    dissolved_per_volume = b[1, c] * rs * rho[2, c] / b[2, c]
    v = vol(1)
    dis = v > 0 ? dissolved_per_volume * v : zero(v)
    brine = v > 0 ? brine - dis : brine - dissolved_per_volume * v
    return (brine, co2 + dis, dis)
end

"""
    ve_boundary_mass_fluxes(bcs, state) -> Matrix{Float64}

Mass fluxes out of the domain (kg/s, positive out) through each boundary
condition in `bcs`, evaluated on `state`, which must hold `:Pressure`,
`:CapillaryPressure`, `:RelativePermeabilities`, `:PhaseViscosities` and
`:PhaseMassDensities` (add them to the reservoir model's
`output_variables` to get them in simulation output); with dissolution and
`dissolved = :carry` also `:ShrinkageFactors` and `:Rs`.

Returns a new 3 × `length(bcs)` matrix: row 1 brine, row 2 CO2 (free and
dissolved), row 3 the dissolved part of row 2. With co2lab's booking of
dissolved CO2 row 1 holds the brine phase's whole mass and row 3 is zero.
Arguments are only read.
"""
function ve_boundary_mass_fluxes(bcs::AbstractVector{<:VEBoundaryCondition}, state)
    q = zeros(3, length(bcs))
    for (k, bc) in enumerate(bcs)
        q_k = _ve_bc_mass_outflow(bc, state, Jutul.gravity_constant)
        q[:, k] .= Jutul.value.(q_k)
    end
    return q
end

# --------------------------------- Jutul hook ---------------------------------

# Modifies: `acc` (the mass-conservation accumulation, per component and cell).
# Component rows: brine, CO2 (free and dissolved).
function Jutul.apply_forces_to_equation!(
        acc, storage, model::CO2LabModel, eq::JutulDarcy.ConservationLaw{:TotalMasses},
        eq_s, force::AbstractVector{<:VEBoundaryCondition}, time)
    state = storage.state
    for bc in force
        q = _ve_bc_mass_outflow(bc, state, Jutul.gravity_constant)
        acc[1, bc.cell] += q[1]
        acc[2, bc.cell] += q[2]
    end
    return acc
end

# Dissolved CO2 carried out by brine, in the dissolved-CO2 law (`:carry`
# only; co2lab's law sees no boundary flux).
# Modifies: `acc` (the dissolved-CO2 residual, per cell).
function Jutul.apply_forces_to_equation!(
        acc, storage, model::CO2LabModel, eq::DissolvedCO2Equation,
        eq_s, force::AbstractVector{<:VEBoundaryCondition}, time)
    state = storage.state
    for bc in force
        bc.carry_dissolved || continue
        acc[1, bc.cell] += _ve_bc_mass_outflow(bc, state, Jutul.gravity_constant)[3]
    end
    return acc
end

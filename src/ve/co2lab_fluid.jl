# co2lab's fluid system: JutulDarcy's black-oil system with co2lab's
# formulation swapped in (roadmap standing decision "Fluid system", Stage 2
# item 0). Moved here from scripts/prototype_blackoil_co2lab.jl.
#
# co2lab's model (CO2VEBlackOilTypeModel) is a black-oil formulation: brine
# plays black-oil's liquid phase, CO2 its vapour phase, and CO2 dissolves in
# brine (Rs). What co2lab does differently is its choice of unknowns and, with
# rate-driven dissolution, one extra conservation law. JutulDarcy's black-oil
# masses, component fluxes and wells read Rs from the state and do not care
# whether it is an unknown, so the system type stays StandardBlackOilSystem,
# re-tagged with formulation :co2lab.
#
# Unknowns follow the configuration, as in MRST (CO2VEBlackOilTypeModel.m
# :19-20, :440-460):
#   always              p, sG
#   dissolution on      + rs, and the dissolved-CO2 conservation law
#   hysteresis on       + sGmax, with co2lab's equation (HysteresisEquation)
#   (res_gas > 0)
# Dissolution off is black-oil without dissolved gas (rs_max = nothing), which
# JutulDarcy already supports in its fluxes and wells.
#
# sGmax always exists. Without hysteresis it is a stored field (a Jutul
# parameter) updated after each converged step: max(sGmax, sG), capped at
# 1 - res_water; it does not affect the flow then (roadmap Stage 2 item 2).
# With hysteresis it is an unknown of the reservoir (roadmap Stage 3).
#
# co2lab's VE relperm, pcWG and PVT are in ve_fluid.jl (Stage 2 items 1, 3).
# Instantaneous dissolution (dis_rate 0 or Inf) is in instant_dissolution.jl.
#
# Note: docstrings cannot be attached to @jutul_secondary definitions; their
# "Modifies:" lines are comments instead.

import JutulDarcy: StandardBlackOilSystem, BlackOilSystem, reference_densities,
                   PhaseVariables, Saturations, TotalMasses,
                   darcy_permeability_potential_differences, phase_upwind

# ----------------------------------- System -----------------------------------

"All systems with co2lab's formulation (two phases, no vaporised brine)."
const CO2LabSystem = StandardBlackOilSystem{<:Any, Nothing, <:Any, <:Any, :co2lab}
"co2lab's formulation with dissolution off (no dissolved gas: `rs_max = nothing`)."
const CO2LabNoDissolutionSystem =
    StandardBlackOilSystem{Nothing, Nothing, <:Any, <:Any, :co2lab}
const CO2LabModel = Jutul.SimulationModel{<:Any, <:CO2LabSystem}
const CO2LabNoDissolutionModel =
    Jutul.SimulationModel{<:Any, <:CO2LabNoDissolutionSystem}

"Saturated dissolved-gas ratio, constant (co2lab's `rsSat` / `dis_max`)."
struct ConstantRs
    value::Float64  # saturated rs
end
(c::ConstantRs)(p) = c.value

"""
    has_dissolution(sys) -> Bool

Whether `sys` carries dissolved CO2 (`rs`).
"""
has_dissolution(sys::StandardBlackOilSystem) = JutulDarcy.has_disgas(sys)

"""
    co2lab_system(; reference_densities = (1020.0, 1.87), rs_sat = nothing,
                    formulation = :co2lab) -> StandardBlackOilSystem

Two-phase black-oil system for co2lab: brine as black-oil's liquid phase, CO2
as its vapour phase. `rs_sat` (surface volume of CO2 per surface volume of
brine) turns dissolution on; `nothing` turns it off.

`formulation = :varswitch` gives the same fluid with JutulDarcy's equilibrium
formulation, used for the wellbores when dissolution is on (see
[`setup_co2lab_model`](@ref)).

Returns a new `StandardBlackOilSystem` with phases `(LiquidPhase(),
VaporPhase())` and the given formulation tag.
"""
function co2lab_system(; reference_densities = (1020.0, 1.87), rs_sat = nothing,
                         formulation::Symbol = :co2lab)
    rs_max = isnothing(rs_sat) ? nothing : ConstantRs(rs_sat)
    # the constructor takes only :varswitch or :zg; other tags are set below
    sys = StandardBlackOilSystem(rs_max = rs_max, rv_max = nothing,
                                 phases = (LiquidPhase(), VaporPhase()),
                                 reference_densities = collect(reference_densities),
                                 formulation = :varswitch, eps_rs = 1e-8)
    return formulation == :varswitch ? sys : _with_formulation(sys, formulation)
end

# StandardBlackOilSystem's constructor only admits :varswitch / :zg
# (JutulDarcy/src/types.jl:174); the struct does not restrict the formulation
# parameter, so rebuild it with another tag.
function _with_formulation(sys::StandardBlackOilSystem, f::Symbol)
    P = typeof(sys).parameters
    T = StandardBlackOilSystem{P[1], P[2], P[3], P[4], f, P[6], P[7], P[8]}
    return T((getfield(sys, k) for k in fieldnames(typeof(sys)))...)
end

# ---------------------------------- Unknowns ----------------------------------

struct GasSaturation <: Jutul.ScalarVariable end
Jutul.minimum_value(::GasSaturation) = 0.0
Jutul.maximum_value(::GasSaturation) = 1.0
Jutul.absolute_increment_limit(::GasSaturation) = 0.2

struct DissolvedGasRatio <: Jutul.ScalarVariable end
Jutul.minimum_value(::DissolvedGasRatio) = 0.0

"""
Jutul hook: co2lab's unknowns, by configuration.

Modifies: `S` (adds the unknowns).
"""
function Jutul.select_primary_variables!(S, sys::CO2LabSystem,
                                         model::Jutul.SimulationModel)
    S[:Pressure] = JutulDarcy.Pressure()
    S[:sG] = GasSaturation()
    if has_dissolution(sys)
        S[:Rs] = DissolvedGasRatio()     # the key black-oil's mass and flux code reads
    end
end

# ----------------------------------- sGmax ------------------------------------

"""
    MaxGasSaturation(res_water, res_gas = 0.0)

co2lab's `sGmax`, between 0 and `1 - res_water`. Without hysteresis a Jutul
parameter, updated explicitly before each step and after each converged
step: `sGmax = min(max(sGmax, sG), 1 - res_water)`. MRST does the same in
`updateState` (`CO2VEBlackOilTypeModel.m:352-353`), but after every Newton
update, so an overshooting iterate can raise its `sGmax`; here only converged
states count. With hysteresis an unknown instead, governed by
[`HysteresisEquation`](@ref), and after every Newton update projected as
co2lab's `updateState` (`:349-357`): `sG ≤ 1 − res_water`, `sGmax` into
`[sG, 1 − res_water]`, then `sG ≥ sGmax·res_gas` (CO2 cannot fall below the
residual its maximum left).
"""
struct MaxGasSaturation <: Jutul.ScalarVariable
    res_water::Float64  # residual brine: sGmax ≤ 1 − res_water
    res_gas::Float64    # residual CO2: sG ≥ sGmax·res_gas (hysteresis)
end
MaxGasSaturation(res_water) = MaxGasSaturation(res_water, 0.0)
Jutul.default_value(model, ::MaxGasSaturation) = 0.0
Jutul.minimum_value(::MaxGasSaturation) = 0.0
Jutul.maximum_value(v::MaxGasSaturation) = 1 - v.res_water
Jutul.absolute_increment_limit(::MaxGasSaturation) = 0.2

"""
Jutul hook (hysteresis: `sGmax` an unknown): the Newton update of `sGmax`,
then co2lab's projection of `sG` and `sGmax` (`updateState`,
`CO2VEBlackOilTypeModel.m:349-357`; see [`MaxGasSaturation`](@ref)). Without
it Newton cycles where a cell switches between drainage and imbibition and an
iterate drops below its residual CO2 (`geothermGradientExample`, a 10-year
step). Where `sG` is not an unknown (instantaneous dissolution) the plain
update.

Modifies: `state[state_symbol]` (`sGmax`) and `state.sG`.
"""
function Jutul.update_primary_variable!(state, v::MaxGasSaturation, state_symbol, model,
                                        dx, w)
    invoke(Jutul.update_primary_variable!,
           Tuple{Any, Jutul.JutulVariables, Any, Any, Any, Any},
           state, v, state_symbol, model, dx, w)
    smax = state[state_symbol]
    haskey(model.primary_variables, :sG) || return smax
    sG, rw, rg = state.sG, v.res_water, v.res_gas
    val = Jutul.value
    for i in eachindex(smax)
        # in co2lab's order: sG capped, sGmax into [sG, 1 − rw], sG above residual
        sg = min(val(sG[i]), 1 - rw)
        m = max(min(val(smax[i]), 1 - rw), sg)
        smax[i] = replace_value(smax[i], m)
        sG[i] = replace_value(sG[i], max(sg, m * rg))
    end
    return smax
end


"""
    _update_max_gas_saturation!(sGmax, v, sG)

Modifies: `sGmax`. `sG` is only read.
"""
function _update_max_gas_saturation!(sGmax, v::MaxGasSaturation, sG)
    for i in eachindex(sGmax, sG)
        # value: sG may carry derivatives, the stored sGmax must not
        sGmax[i] = min(max(sGmax[i], Jutul.value(sG[i])), 1 - v.res_water)
    end
    return sGmax
end

"""
Jutul hook: raise `sGmax` to the saturation the step starts from (this
also covers an initial state with free CO2).

Modifies: `sGmax` (the parameter array).
"""
function Jutul.update_parameter_before_step!(sGmax, v::MaxGasSaturation, storage, model,
                                             dt, forces)
    return _update_max_gas_saturation!(sGmax, v, storage.state.sG)
end

"""
Jutul hook, after a converged step: raise `sGmax` to the converged saturation.

Modifies: `storage.state.sGmax` (the parameter array, shared with
`storage.state0`), for submodels that have it.
"""
function Jutul.update_after_step!(storage, sys::CO2LabSystem, model, dt, forces;
                                  kwarg...)
    prm = storage.variable_definitions.parameters
    # sGmax is a parameter only without hysteresis (else an unknown)
    if haskey(prm, :sGmax)
        _update_max_gas_saturation!(storage.state.sGmax, prm[:sGmax], storage.state.sG)
    end
end

# ----------------------------- Derived quantities -----------------------------

"b = b_ref (1 + c (p − p_ref)) per phase; b = reservoir density / surface density."
struct LinearShrinkage <: PhaseVariables
    b_ref::NTuple{2, Float64}  # b at p_ref, per phase
    c::NTuple{2, Float64}      # compressibility (1/Pa), per phase
    p_ref::NTuple{2, Float64}  # reference pressure (Pa), per phase
end
LinearShrinkage(b_ref, c, p_ref::Real) = LinearShrinkage(b_ref, c, (p_ref, p_ref))
struct BlackOilDensities <: PhaseVariables end
struct ConstantViscosities <: PhaseVariables
    mu::NTuple{2, Float64}  # viscosities (Pa·s), per phase
end
"Mass of CO2 dissolved in brine, per cell."
struct DissolvedCO2Mass <: Jutul.ScalarVariable end

# placeholder PVT until co2lab's is matched
# b_ref: brine 1, CO2 700/1.87 (reservoir/surface density); c per Pa; p_ref 100 bar
_default_shrinkage() = LinearShrinkage((1.0, 700.0 / 1.87), (4.3e-5 / 1e5, 1e-3 / 1e5),
                                       100e5)
_default_viscosity() = ConstantViscosities((8e-4, 6e-5))

"""
Jutul hook: black-oil's derived quantities, minus the two that belong to its
own unknowns, with the placeholder PVT and (with dissolution) the
dissolved-CO2 mass.

Modifies: `S` (the secondary-variable set).
"""
function Jutul.select_secondary_variables!(S, sys::CO2LabSystem,
                                           model::Jutul.SimulationModel)
    # start from black-oil's own set (the generic BlackOilSystem method)
    invoke(Jutul.select_secondary_variables!, Tuple{Any, BlackOilSystem, Any},
           S, sys, model)

    delete!(S, :PhaseState)     # only read by the :varswitch / :zg code
    delete!(S, :Rs)             # an unknown here, when present at all

    S[:ShrinkageFactors]   = _default_shrinkage()
    S[:PhaseMassDensities] = BlackOilDensities()
    S[:PhaseViscosities]   = _default_viscosity()

    # the dissolved-CO2 law, and so its mass, belongs to the reservoir only
    is_well = Jutul.physical_representation(model.domain) isa JutulDarcy.WellDomain
    if has_dissolution(sys) && !is_well
        S[:DissolvedCO2] = DissolvedCO2Mass()
    end
end

# Modifies: `s` (phase saturations).
@jutul_secondary function update_co2lab_saturations!(s, var::Saturations,
                                                    model::CO2LabModel, sG, ix)
    for i in ix
        s[1, i] = 1 - sG[i]
        s[2, i] = sG[i]
    end
end

# Modifies: `b` (shrinkage factors).
@jutul_secondary function update_linear_shrinkage!(b, v::LinearShrinkage, model,
                                                  Pressure, ix)
    for i in ix, ph in 1:2
        b[ph, i] = v.b_ref[ph] * (1 + v.c[ph] * (Pressure[i] - v.p_ref[ph]))
    end
end

# Black-oil densities: brine carries its dissolved CO2.
# Modifies: `rho` (phase mass densities).
@jutul_secondary function update_bo_densities!(rho, var::BlackOilDensities, model,
                                              ShrinkageFactors, Rs, ix)
    ρL, ρV = reference_densities(model.system)
    for i in ix
        rho[1, i] = ShrinkageFactors[1, i] * (ρL + Rs[i] * ρV)
        rho[2, i] = ShrinkageFactors[2, i] * ρV
    end
end

# Densities without dissolution.
# Modifies: `rho` (phase mass densities).
@jutul_secondary function update_bo_densities!(rho, var::BlackOilDensities,
                                              model::CO2LabNoDissolutionModel,
                                              ShrinkageFactors, ix)
    ρS = reference_densities(model.system)
    for i in ix, ph in 1:2
        rho[ph, i] = ShrinkageFactors[ph, i] * ρS[ph]
    end
end

# Modifies: `mu` (phase viscosities).
@jutul_secondary function update_const_visc!(mu, v::ConstantViscosities, model,
                                            Pressure, ix)
    for i in ix, ph in 1:2
        # 0·p: the constant in Pressure's number type (AD or plain)
        mu[ph, i] = v.mu[ph] + 0 * Pressure[i]
    end
end

# Total masses without dissolution. JutulDarcy has methods for "no vaporised
# oil" and for "no dissolved gas" but not for both at once (the two would be
# ambiguous), so the co2lab case gets its own.
# Modifies: `totmass`.
@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses,
                                              model::CO2LabNoDissolutionModel,
                                              ShrinkageFactors, PhaseMassDensities,
                                              Saturations, FluidVolume, ix)
    sys = model.system
    rhoS = reference_densities(sys)
    ind = JutulDarcy.phase_indices(sys)

    for cell in ix
        # Rs, Rv = missing: no dissolved gas, no vaporised oil; false: no water phase
        @views JutulDarcy.blackoil_mass!(totmass[:, cell], FluidVolume,
                                         PhaseMassDensities, missing, missing,
                                         ShrinkageFactors, Saturations, rhoS, cell, ind,
                                         false)
    end
end

# Dissolved part of black-oil's CO2 mass (total_masses.jl, blackoil_mass!):
# Φ · ρ_VS · b_L · s_L · rs.
# Modifies: `m`.
@jutul_secondary function update_dissolved_co2!(m, var::DissolvedCO2Mass,
                                               model::CO2LabModel, Rs, ShrinkageFactors,
                                               Saturations, FluidVolume, ix)
    ρV = reference_densities(model.system)[2]
    for i in ix
        m[i] = FluidVolume[i] * ρV * ShrinkageFactors[1, i] * Saturations[1, i] * Rs[i]
    end
end

# ----------------------------- Dissolved-CO2 law ------------------------------

"""
    DissolvedCO2Equation(flow_discretization, dis_rate, dis_max,
                         vertical = SharpInterface())

Conservation of dissolved CO2 with co2lab's transfer term (rate-driven mode,
`CO2VEBlackOilTypeModel.m:84-97`): `d/dt(M_dis) + div(F_dis) − η = 0`.
Assembled like Jutul's TPFA conservation law, with the same stencil; the flux
is the dissolved part of black-oil's CO2 flux (`JutulDarcy/src/blackoil/
flux.jl`): `ρ_VS · rs_up · (bλ)_L,up · ∇ψ_L`. The transfer `η` is co2lab's
`computeDissolutionTransfer` ([`_dissolution_transfer`](@ref)), with the
brine CO2 has touched from the column-specific vertical model `vertical`.
"""
struct DissolvedCO2Equation{D, M} <: Jutul.JutulEquation
    flow_discretization::D  # faces and stencil (TPFA)
    dis_rate::Float64       # co2lab dis_rate
    dis_max::Float64        # saturated rs
    vertical::M             # vertical model (column-specific): the touched brine
end
DissolvedCO2Equation(flow, dis_rate, dis_max) =
    DissolvedCO2Equation(flow, dis_rate, dis_max, SharpInterface())
Jutul.number_of_equations_per_entity(::Jutul.SimulationModel, ::DissolvedCO2Equation) = 1
function Jutul.local_discretization(eq::DissolvedCO2Equation, i)
    # the half-faces of cell i
    d = eq.flow_discretization
    return view(d.conn_data, d.conn_pos[i]:(d.conn_pos[i+1] - 1))
end

# Declare exactly the conservation law's pattern, in its order
# (Jutul/src/conservation/conservation.jl:487-505): the block matrix layout
# requires every equation's pattern to be identical element by element
# (Jutul/src/models.jl:630), and a traced pattern lists the same entries in a
# different order. The layout also needs the cross-term row-stride fix
# (Jutul.jl #358), which is on the Jutul fork this package uses.
function Jutul.declare_pattern(model, eq::DissolvedCO2Equation, e_s, entity::Jutul.Cells)
    hfd = Array(eq.flow_discretization.conn_data)
    n = Jutul.number_of_entities(model, eq)
    D = collect(1:n)
    I = map(x -> x.self, hfd); J = map(x -> x.other, hfd)

    # neighbour entries first, then each cell's own (the diagonal), as the law
    I, J = Jutul.map_ij_to_active(I, J, model.domain, entity)
    return (vcat(I, D), vcat(J, D))
end

"""
    _co2_face_fluxes(state, model, ldisc) -> (dissolved_out, free_in)

Over the faces of one cell: the net outflow of dissolved CO2 (kg/s, the
brine phase's `ρ_VS · rs · bλ_L · ψ_L`) and the net inflow of free CO2 (the
CO2 phase's `ρ_VS · bλ_V · ψ_V`), both upwinded as black-oil's fluxes.
"""
function _co2_face_fluxes(state, model, ldisc)
    ρV = reference_densities(model.system)[2]
    out_dis = zero(eltype(state.Rs))
    in_free = zero(out_dis)
    for ld in ldisc
        # ψ: phase potential differences over the face, positive out of the cell
        kgrad = Jutul.TPFA(ld.self, ld.other, ld.face_sign)
        upw = Jutul.SPU(ld.self, ld.other)
        ψ = darcy_permeability_potential_differences(ld.face, state, model,
                                                     Jutul.DefaultFlux(), kgrad, upw)
        bλL = phase_upwind(upw, state.SurfaceVolumeMobilities, 1, ψ[1])
        rs = Jutul.upwind(upw, state.Rs, ψ[1])
        out_dis += ρV * rs * bλL * ψ[1]

        # free CO2 counts inflow, hence the sign
        bλV = phase_upwind(upw, state.SurfaceVolumeMobilities, 2, ψ[2])
        in_free -= ρV * bλV * ψ[2]
    end
    return (out_dis, in_free)
end

"""
    _dissolution_transfer(eq, ρV, state, state0, c, dt, out_dis, in_free, src) -> η

co2lab's `computeDissolutionTransfer` (`CO2VEBlackOilTypeModel.m:163-238`) for
cell `c`, as a mass rate (kg/s) — co2lab's `eta .* rhoGS .* bW`. In co2lab's
units (per `ρ_VS · b_W`):

    # co2lab's tie rules (see _max_c in ve_fluid.jl)
    η = _min_v(max(basic_rate, min_transfer), max_transfer)

- `basic_rate = dis_rate · pv / H`;
- `max_transfer`: the lesser of the demand (dissolved CO2 up to `dis_max`
  over the step, plus the dissolved outflow `out_dis`) and the supply (free
  CO2 at the start of the step, plus the free inflow `in_free` and the
  source `src`), at least zero;
- `min_transfer`: brine that has been in contact with CO2 (below the plume,
  above `h_max`) is kept saturated — `dis_max · sW_touched · bW · pv` less
  what was dissolved, at least zero. The untouched brine from the vertical
  model (`_untouched_brine`); for the sharp interface co2lab's.

co2lab's `min_transfer` is in surface volume per time, while the rest is
divided by `bW` once more; kept as co2lab has it.
"""
function _dissolution_transfer(eq, ρV, state, state0, c, dt, out_dis, in_free, src)
    pv, pv0 = state.FluidVolume[c], state0.FluidVolume[c]
    bW, bG0 = state.ShrinkageFactors[1, c], state0.ShrinkageFactors[2, c]
    sW, H = state.Saturations[1, c], state.ColumnHeight[c]

    # most CO2 the brine could still take up, and most CO2 there is to give
    demand = (pv * sW * bW * eq.dis_max * ρV - state0.DissolvedCO2[c]) / dt + out_dis
    supply = pv0 * state0.sG[c] * bG0 * ρV / dt + in_free + src
    max_transfer = _max_c(_min_v(demand, supply), 0.0) / (bW * ρV)

    # brine touched by the plume (above h_max) is kept saturated; sGmax as it
    # stands within the step (as for the heights)
    ρ, smax = state.PhaseMassDensities, max(state.sGmax[c], state.sG[c])
    untouched = _untouched_brine(eq.vertical, c, smax, state.PlumeHeights[2, c], H,
                                 ρ[1, c] - ρ[2, c])
    sW_touched = sW - untouched
    min_dissolved = eq.dis_max * sW_touched * bW * pv
    min_transfer = _max_c(min_dissolved - state0.Rs[c] * sW * bW * pv, 0.0) / dt

    basic_rate = eq.dis_rate * pv / H
    # co2lab's tie rules (see _max_c in ve_fluid.jl)
    η = _min_v(max(basic_rate, min_transfer), max_transfer)
    return η * ρV * bW
end

"""
Jutul hook: the dissolved-CO2 residual of cell `c`.

Modifies: `v`.
"""
function Jutul.update_equation_in_entity!(v, c, state, state0, eq::DissolvedCO2Equation,
                                          model, dt,
                                          ldisc = Jutul.local_discretization(eq, c))
    ρV = reference_densities(model.system)[2]
    acc = (state.DissolvedCO2[c] - state0.DissolvedCO2[c]) / dt
    out_dis, in_free = _co2_face_fluxes(state, model, ldisc)
    # no well source here: DissolutionWellSupplyCT adds its effect
    η = _dissolution_transfer(eq, ρV, state, state0, c, dt, out_dis, in_free, 0.0)
    v[1] = acc + out_dis - η
end

# ---------------------- Well CO2 in the transfer supply -----------------------

"""
    DissolutionWellSupplyCT(reservoir_cells, well_cells)

Cross term from a well to the dissolved-CO2 law: co2lab counts the CO2 a
well puts into its perforated cell, at the current iterate, as supply for
dissolution (`computeDissolutionTransfer`'s `src.value{2}`). In Jutul the
well's inflow reaches the reservoir only through a cross term, so this one
recomputes the transfer of the perforated cell with that inflow (JutulDarcy's
own perforation flux) and adds the difference to the residual, with
derivatives with respect to both reservoir and well. One perforation per
cell, as co2lab's VE wells; several perforations in one cell would each be
counted on their own.
"""
struct DissolutionWellSupplyCT{I <: AbstractVector} <: Jutul.AdditiveCrossTerm
    reservoir_cells::I  # perforated reservoir cells
    well_cells::I       # their well cells
end

Jutul.cross_term_entities(ct::DissolutionWellSupplyCT, eq::DissolvedCO2Equation, model) =
    ct.reservoir_cells
Jutul.cross_term_entities_source(ct::DissolutionWellSupplyCT, eq, model) = ct.well_cells

"""
Jutul hook: the change in the perforated cell's residual from the well's CO2 —
of the dissolved-CO2 law (the transfer) or, with hysteresis, of the
hysteresis equation (the depletion of residual CO2 follows the transfer).

Modifies: `out`.
"""
function Jutul.update_cross_term_in_entity!(out, i, state_t, state0_t,
                                            state_s, state0_s, model_t, model_s,
                                            ct::DissolutionWellSupplyCT, eq, dt,
                                            ldisc = nothing)
    sys = model_t.system
    rhoS = reference_densities(sys)
    conn = JutulDarcy.cross_term_perforation_get_conn(ct, i, state_s, state_t)
    # JutulDarcy's perforation flux, as its ReservoirFromWellFlowCT computes it;
    # untyped buffer: the flux's number type (AD, sparsity tracing) depends on
    # which well quantities it reads (one perforated cell per well, so cheap)
    c = conn.reservoir
    q = Vector{Any}(undef, 2)
    if haskey(state_s, :MassFractions)
        JutulDarcy.simple_well_perforation_flux!(q, sys, state_t, state_s, rhoS, conn)
    else
        JutulDarcy.multisegment_well_perforation_flux!(q, sys, state_t, state_s, rhoS,
                                                       conn)
    end

    # CO2 into the cell (the cross term's convention: positive out of the cell)
    src = -q[2]
    out[1] = _well_supply_change(eq, state_t, state0_t, model_t, c, dt, src)
    return out
end

# Residual with the well's CO2 in the supply, minus without.
function _well_supply_change(eq::DissolvedCO2Equation, state, state0, model, c, dt, src)
    ρV = reference_densities(model.system)[2]
    ld = Jutul.local_discretization(eq, c)
    out_dis, in_free = _co2_face_fluxes(state, model, ld)
    η(q) = _dissolution_transfer(eq, ρV, state, state0, c, dt, out_dis, in_free, q)
    return η(0.0) - η(src)
end


# dissolved-CO2 residual, made dimensionless like CNV: |r| dt / (ρ_VS Φ)
function Jutul.convergence_criterion(model::Jutul.SimulationModel, storage,
                                     eq::DissolvedCO2Equation, eq_s, r;
                                     dt = 1.0, update_report = missing)
    Φ = Jutul.as_value(storage.state.FluidVolume)
    ρV = reference_densities(model.system)[2]
    e = maximum(abs.(vec(r)) .* dt ./ (ρV .* Φ))
    return (CNV = (errors = (e,), names = ("dissolved CO2",)),)
end

# ---------------------------- Hysteresis equation -----------------------------

"""
    HysteresisEquation(flow_discretization, res_water, res_gas, dissolution)

co2lab's equation for `sGmax` with hysteresis on (`res_gas > 0`),
`CO2VEBlackOilTypeModel.m:113-124`, per cell:

    sGmax − max(sGmax0 − s_depletion · dt · fac, sG) = 0

with `sGmax0` from the start of the step. The depletion term, residual CO2
consumed by dissolution, is [`_residual_depletion`](@ref): zero without
dissolution (`dissolution = nothing`), otherwise from the dissolved-CO2 law
`dissolution` ([`DissolvedCO2Equation`](@ref)). `flow_discretization` gives the
matrix pattern (see `declare_pattern` below) and, with dissolution, the faces
the transfer needs.
"""
struct HysteresisEquation{D, E} <: Jutul.JutulEquation
    flow_discretization::D  # faces and matrix pattern
    res_water::Float64      # residual brine saturation
    res_gas::Float64        # residual CO2 saturation
    dissolution::E          # dissolved-CO2 law, or nothing
end
HysteresisEquation(flow, rw, rg) = HysteresisEquation(flow, rw, rg, nothing)
Jutul.number_of_equations_per_entity(::Jutul.SimulationModel, ::HysteresisEquation) = 1
Jutul.local_discretization(::HysteresisEquation{<:Any, Nothing}, i) = nothing
Jutul.local_discretization(eq::HysteresisEquation, i) =
    Jutul.local_discretization(eq.dissolution, i)

# Same pattern as the conservation law, in its order, as for the dissolved-CO2
# law: the block matrix layout requires identical patterns.
function Jutul.declare_pattern(model, eq::HysteresisEquation, e_s, entity::Jutul.Cells)
    law = DissolvedCO2Equation(eq.flow_discretization, 0.0, 0.0)
    return Jutul.declare_pattern(model, law, e_s, entity)
end

"""
    _residual_depletion(eq, state, state0, model, c, dt, ldisc, src) -> Real

The term `s_depletion · dt · fac` of co2lab's hysteresis equation: the drop
in `sGmax` from residual CO2 consumed by dissolution, `fac = (1 − res_water) /
res_gas / (pv · ρ_CO2,ref · b_CO2)`, with `s_depletion` the rate-driven
dissolution transfer (kg/s, [`_dissolution_transfer`](@ref), with the well's
CO2 `src` in its supply) — `CO2VEBlackOilTypeModel.m:95-124`. Zero without
dissolution. (co2lab's instantaneous mode uses the dissolved-CO2 residual
instead: Stage 4 item 2.)
"""
_residual_depletion(eq::HysteresisEquation{<:Any, Nothing}, state, state0, model, c, dt,
                    ldisc, src) = 0.0

function _residual_depletion(eq::HysteresisEquation, state, state0, model, c, dt,
                             ldisc, src)
    ρV = reference_densities(model.system)[2]
    out_dis, in_free = _co2_face_fluxes(state, model, ldisc)
    η = _dissolution_transfer(eq.dissolution, ρV, state, state0, c, dt, out_dis,
                              in_free, src)
    pv, bG = state.FluidVolume[c], state.ShrinkageFactors[2, c]
    fac = (1 - eq.res_water) / eq.res_gas / (pv * ρV * bG)
    return η * dt * fac
end

# The hysteresis residual of cell c, with the well's CO2 src in the transfer.
function _hysteresis_residual(eq, state, state0, model, c, dt, ldisc, src)
    dep = _residual_depletion(eq, state, state0, model, c, dt, ldisc, src)
    return state.sGmax[c] - max(state0.sGmax[c] - dep, state.sG[c])
end

Jutul.cross_term_entities(ct::DissolutionWellSupplyCT, eq::HysteresisEquation, model) =
    ct.reservoir_cells

# Residual with the well's CO2 in the transfer, minus without
# (DissolutionWellSupplyCT).
function _well_supply_change(eq::HysteresisEquation, state, state0, model, c, dt, src)
    ld = Jutul.local_discretization(eq, c)
    return _hysteresis_residual(eq, state, state0, model, c, dt, ld, src) -
           _hysteresis_residual(eq, state, state0, model, c, dt, ld, 0.0)
end

"""
Jutul hook: the hysteresis residual of cell `c` (without well CO2; a
[`DissolutionWellSupplyCT`](@ref) adds it in perforated cells).

Modifies: `v`.
"""
function Jutul.update_equation_in_entity!(v, c, state, state0, eq::HysteresisEquation,
                                          model, dt,
                                          ldisc = Jutul.local_discretization(eq, c))
    v[1] = _hysteresis_residual(eq, state, state0, model, c, dt, ldisc, 0.0)
end

# sGmax residual: a saturation, dimensionless already
function Jutul.convergence_criterion(model::Jutul.SimulationModel, storage,
                                     eq::HysteresisEquation, eq_s, r;
                                     dt = 1.0, update_report = missing)
    return (AbsMax = (errors = (maximum(abs, r),), names = ("sGmax",)),)
end

# ---------------- Black-oil hooks that assume BlackOilUnknown -----------------

function JutulDarcy.handle_alternate_primary_variable_spec!(init, found, rmodel,
                                                            ::CO2LabSystem)
    # co2lab's unknowns are given directly; nothing to translate
    return init
end

function Jutul.convergence_criterion(model::Jutul.SimulationModel{D, S}, storage,
                                     eq::Jutul.ConservationLaw{:TotalMasses}, eq_s, r;
                                     dt = 1.0, update_report = missing
                                     ) where {D, S<:CO2LabSystem}
    # active cells only, without derivatives
    M = Jutul.global_map(model.domain)
    val(x) = Jutul.as_value(Jutul.active_view(x, M, for_variables = false))
    Φ = val(storage.state.FluidVolume)

    b = val(storage.state.ShrinkageFactors)
    rhoS = reference_densities(model.system)
    # black-oil's CNV and MB measures, for two phases
    cnv, mb = JutulDarcy.cnv_mb_errors_bo(r, Φ, b, dt, rhoS, Val(2))
    names = JutulDarcy.phase_names(model.system)
    return (CNV = (errors = cnv, names = names), MB = (errors = mb, names = names))
end

# -------------------------- Model, state, simulator ---------------------------

"""
    setup_co2lab_model(domain; wells = [], fluid = co2lab_ve_fluid(),
                       rock_compressibility = (reference_pressure = 100e5,
                                               compressibility = 1e-10))
        -> MultiModel

VE reservoir model with co2lab's formulation and fluid, plus wells. `domain`
must come from [`ve_domain`](@ref) (it needs `:column_height`; its
`:rugosity` gives the relperm's accretion layer). `fluid`
([`co2lab_ve_fluid`](@ref)) gives the residual saturations, endpoints, PVT
and dissolution. `rock_compressibility` is passed to
[`set_co2lab_rock_compressibility!`](@ref) — co2lab's linear pore-volume
multiplier, with makeVEFluid's defaults; `nothing` keeps the pore volume
constant (e.g. to use JutulDarcy's `set_rock_compressibility!` instead). The
configuration decides the unknowns (see the top of `co2lab_fluid.jl`):

- `fluid.dissolution` turns dissolution on. Rate-driven: `rs` becomes an
  unknown, with co2lab's dissolved-CO2 law ([`DissolvedCO2Equation`](@ref)).
  Instantaneous (`dis_rate` 0 or `Inf`): co2lab's switching unknown `:CO2X`
  ([`DissolvedOrFreeCO2`](@ref)) replaces `sG` and `Rs`, which become
  secondary.
- hysteresis follows co2lab: on when `fluid.res_gas > 0`, making `sGmax` an
  unknown of the reservoir with [`HysteresisEquation`](@ref); with
  dissolution, residual CO2 is consumed by the transfer, as co2lab.
- each well's CO2 counts as supply for dissolution in its cells
  ([`DissolutionWellSupplyCT`](@ref)).
- stacked columns (`ve_columns(...; below_holes = :stack)`, `domain`'s
  `:pillar_columns` above 1), fine 3D cells and discontinuous VE faces
  (`ve_domain(...; fault_flux = :virtual)`) need the transition flux of
  Stage 12: throws.

With dissolution on, the wellbores use JutulDarcy's equilibrium formulation
of the same fluid, while the reservoir uses co2lab's (the dissolved-CO2 law
is not defined on well segments). `setup_reservoir_model`'s own
`wells_systems` argument cannot do this — it is ignored (JutulDarcy.jl #284)
— so the model is assembled here. `domain`, `wells` and `fluid` are only
read.

Returns a new `MultiModel` with submodels `:Reservoir`, one per well (named
as the wells) and `:Facility`. The reservoir submodel has co2lab's VE relperm,
`pcWG`, PVT and (unless `nothing`) pore-volume multiplier, the parameters
`:ColumnHeight` and `:Rugosity`, `:sGmax` (a parameter, or with hysteresis an
unknown with the equation `:hysteresis`), and outputs `:sG`, `:sGmax`,
`:PlumeHeights`, `:CapillaryPressure` (plus `:Rs`, `:DissolvedCO2` with
dissolution) besides JutulDarcy's defaults. The wells get the same PVT.
"""
function setup_co2lab_model(domain; wells = [], fluid::CO2LabVEFluid = co2lab_ve_fluid(),
                            rock_compressibility = (reference_pressure = 100e5,
                                                    compressibility = 1e-10))
    haskey(domain, :pillar_columns) && any(>(1), domain[:pillar_columns]) &&
        error("stacked columns connect through partial faces, which need the " *
              "transition flux of roadmap Stage 12; not implemented yet")
    haskey(domain, :is_ve) && !all(domain[:is_ve]) &&
        error("fine 3D cells need the transition flux of roadmap Stage 12; " *
              "not implemented yet")
    haskey(domain, :face_category) &&
        any(c -> c in _TRANSITIONS, domain[:face_category]) &&
        error("discontinuous VE faces need the transition flux of roadmap " *
              "Stage 12; not implemented yet (fault_flux = :co2lab)")
    # the column-specific vertical model, built once (fringe tables)
    m = _vertical_model(fluid.vertical, domain, fluid)
    model = _co2lab_multimodel(domain, wells, fluid, m)
    rmodel = model[:Reservoir]
    # instantaneous dissolution: co2lab's switching unknown replaces sG and rs
    instantaneous_dissolution(fluid) && _set_instant_dissolution!(rmodel, fluid)

    _set_max_gas_saturation!(rmodel, fluid)
    _add_well_supply_terms!(model)
    push!(rmodel.output_variables, :sG, :sGmax)
    # VE relperm, pcWG, plume heights, column height and rugosity
    _set_ve_fluid!(rmodel, fluid, m)
    isnothing(rock_compressibility) ||
        set_co2lab_rock_compressibility!(rmodel; rock_compressibility...)

    shrinkage, viscosity, well_viscosity = _co2lab_pvt(fluid)
    for m in values(model.models)
        S = m.secondary_variables
        haskey(S, :ShrinkageFactors) || continue
        # the wellbores' equilibrium formulation starts from JutulDarcy's
        # three-phase deck tables; give every submodel the same PVT (wells'
        # viscosities at the pressure: they have no capillary pressure)
        S[:ShrinkageFactors] = shrinkage
        S[:PhaseMassDensities] = BlackOilDensities()
        is_well = JutulDarcy.model_or_domain_is_well(m)
        S[:PhaseViscosities] = is_well ? well_viscosity : viscosity
    end
    needs_temperature(fluid) && _set_temperature!(model, domain)
    return model
end

"""
    _set_temperature!(model, domain)

A `Temperature` parameter for tabulated properties: in the reservoir the
domain's (`ve_domain(...; temperature)`, JutulDarcy's 30 °C where none is
given), in each well its perforated column's (JutulDarcy's wells have none).

Modifies: `model`'s submodels' `parameters` (adds `:Temperature`), the wells'
`data_domain` (sets `:temperature`). `domain` is only read.
"""
function _set_temperature!(model, domain)
    T_res = haskey(domain, :temperature) ? domain[:temperature] :
            fill(Jutul.default_value(nothing, JutulDarcy.Temperature()),
                 number_of_cells(physical_representation(domain)))
    for m in values(model.models)
        haskey(m.secondary_variables, :ShrinkageFactors) || continue
        m.parameters[:Temperature] = JutulDarcy.Temperature()
        d = m.data_domain
        if JutulDarcy.model_or_domain_is_well(m)
            # the wellbore at its (first) perforated column's temperature
            perf = physical_representation(d).perforations.reservoir
            d[:temperature, Cells()] = fill(T_res[first(perf)], number_of_cells(d))
        end
    end
    return model
end

"""
    _set_max_gas_saturation!(rmodel, fluid)

`sGmax` for the reservoir model: a parameter without hysteresis; with
hysteresis (`res_gas > 0`) an unknown, with co2lab's equation.

Modifies: `rmodel.parameters` or `rmodel.primary_variables` and
`rmodel.equations`. `fluid` is only read.
"""
function _set_max_gas_saturation!(rmodel, fluid::CO2LabVEFluid)
    v = MaxGasSaturation(fluid.res_water, fluid.res_gas)
    if fluid.res_gas > 0
        rmodel.primary_variables[:sGmax] = v
        flow = rmodel.domain.discretizations.mass_flow
        # residual CO2 consumed by dissolution: the transfer of the dissolved-CO2
        # law, or with instantaneous dissolution the dissolved-CO2 balance
        dis = instantaneous_dissolution(fluid) ? InstantDissolution() :
              get(rmodel.equations, :dissolved_co2, nothing)
        rmodel.equations[:hysteresis] = HysteresisEquation(flow, fluid.res_water,
                                                           fluid.res_gas, dis)
    else
        rmodel.parameters[:sGmax] = v
    end
    return rmodel
end

"""
    _co2lab_multimodel(domain, wells, fluid, vertical) -> MultiModel

The reservoir (co2lab's formulation; with dissolution also the dissolved-CO2
law, with the column-specific vertical model `vertical`) and the wells (with
dissolution on JutulDarcy's equilibrium formulation). Arguments are only
read.
"""
function _co2lab_multimodel(domain, wells, fluid::CO2LabVEFluid, vertical)
    rho_ref = fluid.rho_ref
    # nothing: no dissolved CO2 at all (rs is then not even a variable)
    rs_sat = fluid.dissolution ? fluid.dis_max : nothing
    sys = co2lab_system(reference_densities = rho_ref, rs_sat = rs_sat)
    # without dissolution one system serves reservoir and wells alike
    isnothing(rs_sat) && return setup_reservoir_model(domain, sys; wells = wells)

    # with it, the wellbores need JutulDarcy's own (variable-switching) formulation
    sys_well = co2lab_system(reference_densities = rho_ref, rs_sat = rs_sat,
                             formulation = :varswitch)
    full = setup_reservoir_model(domain, sys_well; wells = wells)

    # the reservoir alone, in co2lab's formulation, plus the dissolved-CO2 law
    rmodel = setup_reservoir_model(domain, sys)[:Reservoir]
    flow = rmodel.domain.discretizations.mass_flow
    rmodel.equations[:dissolved_co2] = DissolvedCO2Equation(flow, fluid.dis_rate,
                                                            fluid.dis_max, vertical)
    push!(rmodel.output_variables, :Rs, :DissolvedCO2)

    # the wells and facility from `full`, the reservoir with co2lab's formulation
    models = (; (k => (k == :Reservoir ? rmodel : full.models[k])
                 for k in keys(full.models))...)
    return MultiModel(models; cross_terms = full.cross_terms, groups = full.groups,
                      context = full.context, reduction = full.reduction)
end

"""
    _add_well_supply_terms!(model)

Each well's CO2 counts as supply for dissolution in its perforated cells
([`DissolutionWellSupplyCT`](@ref)): one cross term per well for the
dissolved-CO2 law and, with hysteresis, one for the hysteresis equation.

Modifies: `model` (its cross terms).
"""
function _add_well_supply_terms!(model)
    eqs = model[:Reservoir].equations
    for (k, m) in pairs(model.models), eq in (:dissolved_co2, :hysteresis)
        g = Jutul.physical_representation(m.domain)
        # wells only, and only where a rate-driven transfer exists
        g isa JutulDarcy.WellDomain && haskey(eqs, eq) && haskey(eqs, :dissolved_co2) ||
            continue
        perf = g.perforations
        ct = DissolutionWellSupplyCT(vec(perf.reservoir), vec(perf.self))
        Jutul.add_cross_term!(model, ct, target = :Reservoir, source = k, equation = eq)
    end
    return model
end

"""
    setup_co2lab_state(model; pressure, sG = 0.0, rs = 0.0, forces = nothing) -> state

Initial state for a model from [`setup_co2lab_model`](@ref). `pressure` is
per reservoir cell (or a scalar). `rs` is only used with dissolution on;
with hysteresis `sGmax` starts at `sG`, as co2lab's examples. `model` and
`forces` are only read.

Well nodes start at their first perforated cell's pressure, filled with what
the well injects according to `forces` (from `setup_reservoir_forces`, or a
vector of them: the first is used) — CO2 for `co2_injector_control` — and
with brine otherwise. co2lab's wells have no wellbore volume; a wellbore
starting with brine would push it into the reservoir with the first CO2
(measured on the sloping aquifer: 1e-5 in `sG` and BHP, 1e-4 in pressure
with instantaneous dissolution; with the wellbore filled with the injected
CO2, the well runs agree with co2lab to solver precision).

Written per submodel: `setup_reservoir_state` copies the reservoir's unknowns
to the wells (`JutulDarcy/src/state0.jl:297-300`), which does not hold when
the wellbores use another formulation.

Returns a new state as `Jutul.setup_state` builds it: a `Dict` with one entry
per submodel (`:Reservoir`, the wells, `:Facility`), each a `Dict{Symbol, Any}`
of that submodel's unknowns.
"""
function setup_co2lab_state(model; pressure, sG = 0.0, rs = 0.0, forces = nothing)
    rmodel = model[:Reservoir]
    nc = Jutul.number_of_cells(rmodel.domain)
    p = pressure isa Number ? fill(Float64(pressure), nc) : collect(Float64, pressure)
    init = Dict{Symbol, Any}(:Reservoir => _reservoir_init(rmodel, p, sG, rs, nc))

    perforated_cells(wm) =
        Jutul.physical_representation(wm.domain).perforations.reservoir
    first_perf_pressure(wm) = p[first(perforated_cells(wm))]

    for (k, m) in pairs(model.models)
        k == :Reservoir && continue
        g = Jutul.physical_representation(m.domain)
        if g isa JutulDarcy.WellGroup
            # the facility: no rates yet, each BHP at its well's first perforation
            p_bh = [first_perf_pressure(model[w]) for w in g.well_symbols]
            init[k] = Dict(:TotalSurfaceMassRate => 0.0, :SurfacePhaseRates => 0.0,
                           :BottomHolePressure => p_bh)
        else
            init[k] = _wellbore_init(m, first_perf_pressure(m), _injected_sg(forces, k))
        end
    end
    return setup_state(model, init)
end

# The reservoir's unknowns: pressure and sG, plus rs (dissolution), sGmax
# (hysteresis) or co2lab's switching unknown (instantaneous dissolution).
function _reservoir_init(rmodel, p, sG, rs, nc)
    P = rmodel.primary_variables
    rinit = Dict{Symbol, Any}(:Pressure => p, :sG => sG)
    has_dissolution(rmodel.system) && (rinit[:Rs] = rs)
    haskey(P, :sGmax) && (rinit[:sGmax] = sG)
    # instantaneous dissolution: one unknown X per cell instead of sG and rs
    if haskey(P, :CO2X)
        v = P[:CO2X]
        sGv, rsv = (x isa Number ? fill(Float64(x), nc) : x for x in (sG, rs))
        rinit[:CO2X] = [_co2lab_x(v, sGv[i], rsv[i]) for i in 1:nc]
        delete!(rinit, :sG); delete!(rinit, :Rs)
    end
    return rinit
end

# CO2 saturation of what well `k` injects under `forces` (0 if it does not
# inject, or no forces are given): the CO2 phase's share of its control's
# wellbore mixture (InjectorControl's `phases`).
function _injected_sg(forces, k)
    isnothing(forces) && return 0.0
    f = forces isa AbstractVector ? first(forces) : forces
    ctrl = get(f[:Facility].control, k, nothing)
    ctrl isa JutulDarcy.InjectorControl || return 0.0
    # phases: (phase, fraction) pairs; phase 2 is CO2
    return sum((x[2] for x in ctrl.phases if x[1] == 2); init = 0.0)
end

# A wellbore at pressure pw with CO2 saturation sg (the rest brine, without
# dissolved CO2 unless saturated with free CO2 present).
function _wellbore_init(m, pw, sg)
    if haskey(m.primary_variables, :BlackOilUnknown)
        # equilibrium wellbore: black-oil's switching unknown
        x = sg > 0 ? JutulDarcy.BlackOilX(sg, JutulDarcy.OilAndGas, false) :
                     JutulDarcy.BlackOilX(m.system, pw; sg = 0.0, rs = 0.0)
        return Dict(:Pressure => pw, :BlackOilUnknown => x)
    end
    # co2lab-formulation wellbore (no dissolution): sG itself
    return Dict(:Pressure => pw, :sG => sg)
end

"""
    setup_co2lab_simulator(case; dissolved_co2_tolerance = 1e-6,
                           hysteresis_tolerance = 1e-6, kwarg...) -> (sim, cfg)

`setup_reservoir_simulator`, plus tolerances for co2lab's extra equations
when present: the dissolved-CO2 law (dimensionless, like CNV) and the
hysteresis equation (a saturation; co2lab's default `nonlinearTolerance` is
1e-6 on it). `case` is only read.

Returns `(sim, cfg)`: a new `Simulator` and its `JutulConfig`, with
`cfg[:tolerances][:Reservoir][:dissolved_co2]` and `[:hysteresis]` set when the
model has those equations.
"""
function setup_co2lab_simulator(case; dissolved_co2_tolerance = 1e-6,
                                hysteresis_tolerance = 1e-6, kwarg...)
    sim, cfg = setup_reservoir_simulator(case; kwarg...)

    eqs, tol = case.model[:Reservoir].equations, cfg[:tolerances][:Reservoir]
    haskey(eqs, :dissolved_co2) && (tol[:dissolved_co2] = dissolved_co2_tolerance)
    haskey(eqs, :hysteresis) && (tol[:hysteresis] = hysteresis_tolerance)
    return (sim, cfg)
end

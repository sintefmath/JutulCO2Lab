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
#   hysteresis on       + sGmax as an unknown with its own residual (Stage 3,
#                         not implemented yet)
# Dissolution off is black-oil without dissolved gas (rs_max = nothing), which
# JutulDarcy already supports in its fluxes and wells.
#
# sGmax always exists, as a stored field (a Jutul parameter) updated
# explicitly: max(sGmax, sG), capped at 1 - res_water. Without hysteresis it
# does not affect the flow (roadmap Stage 2 item 2); the trapping inventory
# needs it later.
#
# co2lab's VE relperm, pcWG and PVT are in ve_fluid.jl (Stage 2 items 1, 3).
# Still placeholders: the transfer rate (linear relaxation with a smoothed
# onset, not co2lab's computeDissolutionTransfer — Stage 5), and
# instantaneous dissolution (Stage 5).
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
    value::Float64
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
    MaxGasSaturation(res_water)

co2lab's `sGmax`, stored as a Jutul parameter and updated explicitly before
each step and after each converged step: `sGmax = min(max(sGmax, sG),
1 - res_water)`. MRST does the same in `updateState`
(`CO2VEBlackOilTypeModel.m:352-353`), but after every Newton update, so an
overshooting iterate can raise its `sGmax`; here only converged states count.
With hysteresis (Stage 3) it becomes an unknown instead.
"""
struct MaxGasSaturation <: Jutul.ScalarVariable
    res_water::Float64
end
Jutul.default_value(model, ::MaxGasSaturation) = 0.0

"""
    _update_max_gas_saturation!(sGmax, v, sG)

Modifies: `sGmax`. `sG` is only read.
"""
function _update_max_gas_saturation!(sGmax, v::MaxGasSaturation, sG)
    for i in eachindex(sGmax, sG)
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
    if haskey(prm, :sGmax)
        _update_max_gas_saturation!(storage.state.sGmax, prm[:sGmax], storage.state.sG)
    end
end

# ----------------------------- Derived quantities -----------------------------

"b = b_ref (1 + c (p − p_ref)) per phase; b = reservoir density / surface density."
struct LinearShrinkage <: PhaseVariables
    b_ref::NTuple{2, Float64}
    c::NTuple{2, Float64}
    p_ref::NTuple{2, Float64}
end
LinearShrinkage(b_ref, c, p_ref::Real) = LinearShrinkage(b_ref, c, (p_ref, p_ref))
struct BlackOilDensities <: PhaseVariables end
struct ConstantViscosities <: PhaseVariables
    mu::NTuple{2, Float64}
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
    DissolvedCO2Equation(flow_discretization, rate)

Conservation of dissolved CO2 with a transfer term (co2lab's rate-driven
mode, `CO2VEBlackOilTypeModel.m:84-97`): `d/dt(M_dis) + div(F_dis) − η = 0`.
Assembled like Jutul's TPFA conservation law, with the same stencil; the flux
is the dissolved part of black-oil's CO2 flux (`JutulDarcy/src/blackoil/
flux.jl`): `ρ_VS · rs_up · (bλ)_L,up · ∇ψ_L`.

`rate` (1/s) is a placeholder transfer law: linear relaxation towards
saturation where free CO2 is present, with the onset smoothed over
`sG ∈ [0, 1e-3]` (a hard on/off test made Newton cycle). co2lab's
`computeDissolutionTransfer` replaces it (Stage 5).
"""
struct DissolvedCO2Equation{D} <: Jutul.JutulEquation
    flow_discretization::D
    rate::Float64
end
Jutul.number_of_equations_per_entity(::Jutul.SimulationModel, ::DissolvedCO2Equation) = 1
function Jutul.local_discretization(eq::DissolvedCO2Equation, i)
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

    I, J = Jutul.map_ij_to_active(I, J, model.domain, entity)
    return (vcat(I, D), vcat(J, D))
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

    div = zero(acc)
    for ld in ldisc
        kgrad = Jutul.TPFA(ld.self, ld.other, ld.face_sign)
        upw = Jutul.SPU(ld.self, ld.other)
        ψ = darcy_permeability_potential_differences(ld.face, state, model,
                                                     Jutul.DefaultFlux(), kgrad, upw)
        ψL = ψ[1]
        bλL = phase_upwind(upw, state.SurfaceVolumeMobilities, 1, ψL)
        rs = Jutul.upwind(upw, state.Rs, ψL)
        div += ρV * rs * bλL * ψL
    end

    rs_sat = first(model.system.rs_max)(state.Pressure[c])   # single PVT region
    free = min(max(state.sG[c], 0.0) / 1e-3, 1.0)
    # transfer: rate × rs deficit × CO2 density × brine surface volume (b_L s_L Φ)
    η = free * eq.rate * max(rs_sat - state.Rs[c], 0.0) * ρV *
        state.ShrinkageFactors[1, c] * state.Saturations[1, c] * state.FluidVolume[c]
    v[1] = acc + div - η
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

# ---------------- Black-oil hooks that assume BlackOilUnknown -----------------

function JutulDarcy.handle_alternate_primary_variable_spec!(init, found, rmodel,
                                                            ::CO2LabSystem)
    return init
end

function Jutul.convergence_criterion(model::Jutul.SimulationModel{D, S}, storage,
                                     eq::Jutul.ConservationLaw{:TotalMasses}, eq_s, r;
                                     dt = 1.0, update_report = missing
                                     ) where {D, S<:CO2LabSystem}
    M = Jutul.global_map(model.domain)
    val(x) = Jutul.as_value(Jutul.active_view(x, M, for_variables = false))
    Φ = val(storage.state.FluidVolume)

    b = val(storage.state.ShrinkageFactors)
    rhoS = reference_densities(model.system)
    cnv, mb = JutulDarcy.cnv_mb_errors_bo(r, Φ, b, dt, rhoS, Val(2))
    names = JutulDarcy.phase_names(model.system)
    return (CNV = (errors = cnv, names = names), MB = (errors = mb, names = names))
end

# -------------------------- Model, state, simulator ---------------------------

"""
    setup_co2lab_model(domain; wells = [], fluid = co2lab_ve_fluid(),
                       rs_sat = nothing, transfer_rate = nothing) -> MultiModel

VE reservoir model with co2lab's formulation and fluid, plus wells. `domain`
must come from [`ve_domain`](@ref) (it needs `:column_height`). `fluid`
([`co2lab_ve_fluid`](@ref)) gives the residual saturations, endpoints,
rugosity, PVT and pore-volume multiplier. The configuration decides the
unknowns (see the top of `co2lab_fluid.jl`):

- `rs_sat` and `transfer_rate` (1/s) turn rate-driven dissolution on. Both or
  neither must be given. Instantaneous dissolution is Stage 5 work and throws.
- hysteresis follows co2lab: on when `fluid.res_gas > 0` — Stage 3 work,
  throws for now.

With dissolution on, the wellbores use JutulDarcy's equilibrium formulation
of the same fluid, while the reservoir uses co2lab's (the dissolved-CO2 law
is not defined on well segments). `setup_reservoir_model`'s own
`wells_systems` argument cannot do this — it is ignored (JutulDarcy.jl #284)
— so the model is assembled here. `domain`, `wells` and `fluid` are only
read.

Returns a new `MultiModel` with submodels `:Reservoir`, one per well (named
as the wells) and `:Facility`. The reservoir submodel has co2lab's VE relperm,
`pcWG`, PVT and pore-volume multiplier, the parameters `:sGmax` and
`:ColumnHeight`, and outputs `:sG`, `:sGmax`, `:PlumeHeights`,
`:CapillaryPressure` (plus `:Rs`, `:DissolvedCO2` with dissolution) besides
JutulDarcy's defaults. The wells get the same PVT.
"""
function setup_co2lab_model(domain; wells = [], fluid::CO2LabVEFluid = co2lab_ve_fluid(),
                            rs_sat = nothing, transfer_rate = nothing)
    fluid.res_gas > 0 && error("hysteresis (res_gas > 0: sGmax as an unknown) is " *
                               "roadmap Stage 3; not implemented yet")
    isnothing(rs_sat) == isnothing(transfer_rate) ||
        error("rate-driven dissolution needs both `rs_sat` and `transfer_rate`")
    if !isnothing(rs_sat) && !(0 < transfer_rate < Inf)
        error("instantaneous dissolution (co2lab's dis_rate = 0 or Inf) is roadmap " *
              "Stage 5; not implemented yet")
    end

    model = _co2lab_multimodel(domain, wells, fluid.rho_ref, rs_sat, transfer_rate)
    rmodel = model[:Reservoir]

    rmodel.parameters[:sGmax] = MaxGasSaturation(fluid.res_water)
    push!(rmodel.output_variables, :sG, :sGmax)
    _set_ve_fluid!(rmodel, fluid)

    shrinkage, viscosity = _co2lab_pvt(fluid)
    for m in values(model.models)
        S = m.secondary_variables
        haskey(S, :ShrinkageFactors) || continue
        # the wellbores' equilibrium formulation starts from JutulDarcy's
        # three-phase deck tables; give every submodel the same PVT
        S[:ShrinkageFactors] = shrinkage
        S[:PhaseMassDensities] = BlackOilDensities()
        S[:PhaseViscosities] = viscosity
    end
    return model
end

"""
    _co2lab_multimodel(domain, wells, rho_ref, rs_sat, transfer_rate) -> MultiModel

The reservoir (co2lab's formulation; with dissolution also the dissolved-CO2
law) and the wells (with dissolution on JutulDarcy's equilibrium
formulation). Arguments are only read.
"""
function _co2lab_multimodel(domain, wells, rho_ref, rs_sat, transfer_rate)
    sys = co2lab_system(reference_densities = rho_ref, rs_sat = rs_sat)
    isnothing(rs_sat) && return setup_reservoir_model(domain, sys; wells = wells)

    sys_well = co2lab_system(reference_densities = rho_ref, rs_sat = rs_sat,
                             formulation = :varswitch)
    full = setup_reservoir_model(domain, sys_well; wells = wells)

    rmodel = setup_reservoir_model(domain, sys)[:Reservoir]
    flow = rmodel.domain.discretizations.mass_flow
    rmodel.equations[:dissolved_co2] = DissolvedCO2Equation(flow, transfer_rate)
    push!(rmodel.output_variables, :Rs, :DissolvedCO2)

    # the wells and facility from `full`, the reservoir with co2lab's formulation
    models = (; (k => (k == :Reservoir ? rmodel : full.models[k])
                 for k in keys(full.models))...)
    return MultiModel(models; cross_terms = full.cross_terms, groups = full.groups,
                      context = full.context, reduction = full.reduction)
end

"""
    setup_co2lab_state(model; pressure, sG = 0.0, rs = 0.0) -> state

Initial state for a model from [`setup_co2lab_model`](@ref). `pressure` is
per reservoir cell (or a scalar); well nodes start at their first perforated
cell's pressure, with no gas. `rs` is only used with dissolution on. `model`
is only read.

Written per submodel: `setup_reservoir_state` copies the reservoir's unknowns
to the wells (`JutulDarcy/src/state0.jl:297-300`), which does not hold when
the wellbores use another formulation.

Returns a new state as `Jutul.setup_state` builds it: a `Dict` with one entry
per submodel (`:Reservoir`, the wells, `:Facility`), each a `Dict{Symbol, Any}`
of that submodel's unknowns.
"""
function setup_co2lab_state(model; pressure, sG = 0.0, rs = 0.0)
    rmodel = model[:Reservoir]
    nc = Jutul.number_of_cells(rmodel.domain)
    p = pressure isa Number ? fill(Float64(pressure), nc) : collect(Float64, pressure)
    rinit = Dict{Symbol, Any}(:Pressure => p, :sG => sG)
    has_dissolution(rmodel.system) && (rinit[:Rs] = rs)
    init = Dict{Symbol, Any}(:Reservoir => rinit)

    perforated_cells(wm) =
        Jutul.physical_representation(wm.domain).perforations.reservoir
    first_perf_pressure(wm) = p[first(perforated_cells(wm))]

    for (k, m) in pairs(model.models)
        k == :Reservoir && continue
        g = Jutul.physical_representation(m.domain)
        if g isa JutulDarcy.WellGroup
            p_bh = [first_perf_pressure(model[w]) for w in g.well_symbols]
            init[k] = Dict(:TotalSurfaceMassRate => 0.0, :SurfacePhaseRates => 0.0,
                           :BottomHolePressure => p_bh)
        else
            pw = first_perf_pressure(m)
            if haskey(m.primary_variables, :BlackOilUnknown)
                # equilibrium wellbore: black-oil's switching unknown, no gas
                x = JutulDarcy.BlackOilX(m.system, pw; sg = 0.0, rs = 0.0)
                init[k] = Dict(:Pressure => pw, :BlackOilUnknown => x)
            else
                init[k] = Dict(:Pressure => pw, :sG => 0.0)
            end
        end
    end
    return setup_state(model, init)
end

"""
    setup_co2lab_simulator(case; dissolved_co2_tolerance = 1e-6, kwarg...) -> (sim, cfg)

`setup_reservoir_simulator`, plus a tolerance for the dissolved-CO2
equation (dimensionless, like CNV) when present. `case` is only read.

Returns `(sim, cfg)`: a new `Simulator` and its `JutulConfig`, with
`cfg[:tolerances][:Reservoir][:dissolved_co2]` set when the model has that
equation.
"""
function setup_co2lab_simulator(case; dissolved_co2_tolerance = 1e-6, kwarg...)
    sim, cfg = setup_reservoir_simulator(case; kwarg...)

    if haskey(case.model[:Reservoir].equations, :dissolved_co2)
        cfg[:tolerances][:Reservoir][:dissolved_co2] = dissolved_co2_tolerance
    end
    return (sim, cfg)
end

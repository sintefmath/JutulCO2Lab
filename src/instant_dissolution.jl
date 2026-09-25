# co2lab's instantaneous dissolution (makeVEFluid's dis_rate 0 or Inf; roadmap
# Stage 4 item 2).
#
# Brine is saturated (rs = dis_max) wherever free CO2 is present. co2lab solves
# for one unknown X per cell besides pressure (CO2VEBlackOilTypeModel.m:
# getPrimaryVariables, initStateAD, updateState): X is sG where the brine is
# saturated and rs where it is not (then sG = 0). The switch happens in the
# update after each Newton iteration, classified on the state before it:
#
#   undersaturated: rs += dX, capped to [0, dis_max]; at dis_max the cell is
#                   saturated from then on, with sG = 0
#   saturated:      sG += dX, the step limited to dsMaxAbs = 0.2 and sG capped
#                   to [0, 1 − res_water]; a cell whose sG falls to 0 turns
#                   undersaturated, with rs = dis_max − eps
#
# Here X is JutulDarcy's black-oil switching value BlackOilX (value plus phase
# state: OilOnly = undersaturated, OilAndGas = saturated), under a variable
# type of our own with co2lab's update rule instead of JutulDarcy's (roadmap
# standing decision). sG and Rs become secondary variables computed from X, so
# everything that reads them is unchanged. There is no dissolved-CO2 law; with
# hysteresis the residual CO2 consumed by dissolution is the dissolved-CO2
# balance's residual (co2lab's s_depletion = max(eq_dis, 0)).

import JutulDarcy: BlackOilX, OilOnly, OilAndGas
import Jutul: replace_value

# -------------------------------- The unknown ---------------------------------

"""
    DissolvedOrFreeCO2(dis_max, res_water; ds_max = 0.2)

co2lab's unknown `X` for instantaneous dissolution: `rs` in undersaturated
cells, `sG` in saturated ones, stored as JutulDarcy's `BlackOilX` (phase state
`OilOnly` / `OilAndGas`). Updated with co2lab's switching rule (see the top of
`instant_dissolution.jl`); `ds_max` is co2lab's `dsMaxAbs`.
"""
struct DissolvedOrFreeCO2 <: Jutul.ScalarVariable
    dis_max::Float64
    res_water::Float64
    ds_max::Float64
end
DissolvedOrFreeCO2(dis_max, res_water; ds_max = 0.2) =
    DissolvedOrFreeCO2(Float64(dis_max), Float64(res_water), Float64(ds_max))

"""
    _co2lab_x(v::DissolvedOrFreeCO2, sG, rs) -> BlackOilX

The unknown for a cell with CO2 saturation `sG` and dissolved ratio `rs`:
saturated (`sG`) if there is free CO2, otherwise undersaturated (`rs`).
"""
_co2lab_x(v::DissolvedOrFreeCO2, sG, rs) =
    sG > 0 ? BlackOilX(Float64(sG), OilAndGas, false) :
             BlackOilX(Float64(min(rs, v.dis_max - eps())), OilOnly, false)

# JutulDarcy's storage hooks for BlackOilX, for this variable type
# (JutulDarcy/src/blackoil/variables/variables.jl:105-195)
Jutul.scalarized_primary_variable_type(model, ::DissolvedOrFreeCO2, T = Float64) =
    BlackOilX{T}

function Jutul.scalarize_primary_variable(model, source_vec, ::DissolvedOrFreeCO2, index;
                                          numeric::Bool = false)
    x = source_vec[index]
    xval = Jutul.value(x.val)
    return numeric ? xval : BlackOilX(xval, x.phases_present, x.sat_close)
end

function Jutul.descalarize_primary_variable!(dest_array, model, V::BlackOilX,
                                             ::DissolvedOrFreeCO2, index, ref = missing;
                                             F = identity)
    V_old = dest_array[index]
    val = replace_value(V_old.val, F(V.val))
    dest_array[index] = BlackOilX(val, V.phases_present, V.sat_close)
end

function Jutul.descalarize_primary_variable!(dest_array, model, V, ::DissolvedOrFreeCO2,
                                             index, ref = missing; F = identity)
    V_old = dest_array[index]
    val = replace_value(V_old.val, F(V))
    dest_array[index] = BlackOilX(val, V_old.phases_present, V_old.sat_close)
end

function Jutul.initialize_primary_variable_ad!(state, model, ::DissolvedOrFreeCO2, symb,
                                               npartials; offset, kwarg...)
    pre = state[symb]
    vals = map(x -> x.val, pre)
    ad = Jutul.allocate_array_ad(vals, diag_pos = offset + 1, context = model.context,
                                 npartials = npartials; kwarg...)
    state[symb] = map((v, x) -> BlackOilX(v, x.phases_present, false), ad, pre)
    return state
end

function Jutul.increment_norm(dX, state, model, X, ::DissolvedOrFreeCO2)
    sg = rs = zero(eltype(dX))
    sg_sum = rs_sum = zero(eltype(dX))
    for i in 1:length(dX)
        d = abs(dX[i])
        if X[i].phases_present == OilAndGas
            sg, sg_sum = max(sg, d), sg_sum + d
        else
            rs, rs_sum = max(rs, d), rs_sum + d
        end
    end
    return (sg_sum = sg_sum, sg_max = sg, rs_sum = rs_sum, rs_max = rs)
end

"""
Jutul hook: co2lab's update of `X` after a Newton iteration, with the switch
between undersaturated (`rs`) and saturated (`sG`) cells.

Modifies: `state[state_symbol]` (the unknown).
"""
function Jutul.update_primary_variable!(state, v::DissolvedOrFreeCO2, state_symbol,
                                        model, Dx, w)
    X = state[state_symbol]
    active = Jutul.active_entities(model.domain, Cells(), for_variables = true)
    for (i, dx) in zip(active, Dx)
        X[i] = _update_co2lab_x(v, X[i], w * dx)
    end
    return X
end

# One cell of co2lab's update (CO2VEBlackOilTypeModel.m:403-430).
function _update_co2lab_x(v, x::BlackOilX, dx)
    val = x.val
    if x.phases_present == OilOnly
        rs = Jutul.value(val) + dx
        # at dis_max the brine is saturated; free CO2 starts at zero
        rs >= v.dis_max && return BlackOilX(replace_value(val, 0.0), OilAndGas, false)
        return BlackOilX(replace_value(val, max(rs, 0.0)), OilOnly, false)
    end

    ds = clamp(dx, -v.ds_max, v.ds_max)
    sg = clamp(Jutul.value(val) + ds, 0.0, 1 - v.res_water)
    # free CO2 gone: undersaturated, just below dis_max
    ds < 0 && sg == 0 &&
        return BlackOilX(replace_value(val, v.dis_max - eps()), OilOnly, false)
    return BlackOilX(replace_value(val, sg), OilAndGas, false)
end

# ------------------------------ sG and Rs from X ------------------------------

"`sG` from co2lab's unknown `X`: `X` where saturated, zero elsewhere."
struct InstantGasSaturation <: Jutul.ScalarVariable end

"`rs` from co2lab's unknown `X`: `dis_max` where saturated, `X` elsewhere."
struct InstantDissolvedGasRatio <: Jutul.ScalarVariable
    dis_max::Float64
end

# Modifies: `sG`.
@jutul_secondary function update_instant_sg!(sG, v::InstantGasSaturation, model,
                                            CO2X, ix)
    for i in ix
        x = CO2X[i]
        sG[i] = x.phases_present == OilAndGas ? x.val : zero(x.val)
    end
end

# Modifies: `rs`.
@jutul_secondary function update_instant_rs!(rs, v::InstantDissolvedGasRatio, model,
                                            CO2X, ix)
    for i in ix
        x = CO2X[i]
        rs[i] = x.phases_present == OilAndGas ? v.dis_max + zero(x.val) : x.val
    end
end

# ------------------------------- Model settings -------------------------------

"""
    _set_instant_dissolution!(rmodel, fluid)

Switch a reservoir model built for rate-driven dissolution to co2lab's
instantaneous mode: unknown `:CO2X` ([`DissolvedOrFreeCO2`](@ref)) instead of
`sG` and `Rs` (which become secondary), no dissolved-CO2 law.

Modifies: `rmodel.primary_variables`, `rmodel.secondary_variables`,
`rmodel.equations`. `fluid` is only read.
"""
function _set_instant_dissolution!(rmodel, fluid::CO2LabVEFluid)
    P, S = rmodel.primary_variables, rmodel.secondary_variables
    delete!(P, :sG); delete!(P, :Rs)
    delete!(rmodel.equations, :dissolved_co2)

    # keep pressure first, then X, then whatever else (sGmax with hysteresis)
    rest = [k => P[k] for k in keys(P) if k != :Pressure]
    foreach(k -> delete!(P, k), first.(rest))
    P[:CO2X] = DissolvedOrFreeCO2(fluid.dis_max, fluid.res_water)
    foreach(kv -> P[kv[1]] = kv[2], rest)

    S[:sG] = InstantGasSaturation()
    S[:Rs] = InstantDissolvedGasRatio(fluid.dis_max)
    return rmodel
end

# -------------------------- Hysteresis: residual CO2 --------------------------

"Marker for [`HysteresisEquation`](@ref): instantaneous dissolution."
struct InstantDissolution end

# co2lab's s_depletion for instantaneous dissolution: the dissolved-CO2
# balance's residual where positive (CO2VEBlackOilTypeModel.m:100-107) — the
# dissolved mass gained beyond what the brine brought in.
function _residual_depletion(eq::HysteresisEquation{<:Any, InstantDissolution}, state,
                             state0, model, c, dt, ldisc, src)
    ρV = reference_densities(model.system)[2]
    out_dis, _ = _co2_face_fluxes(state, model, ldisc)
    eq_dis = (state.DissolvedCO2[c] - state0.DissolvedCO2[c]) / dt + out_dis
    pv, bG = state.FluidVolume[c], state.ShrinkageFactors[2, c]
    fac = (1 - eq.res_water) / eq.res_gas / (pv * ρV * bG)
    return _max_c(eq_dis, 0.0) * dt * fac
end

Jutul.local_discretization(eq::HysteresisEquation{<:Any, InstantDissolution}, i) =
    Jutul.local_discretization(DissolvedCO2Equation(eq.flow_discretization, 0.0, 0.0), i)

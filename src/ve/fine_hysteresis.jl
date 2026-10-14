# co2lab's endpoint-scaling hysteresis for fine 3D cells, as JutulDarcy
# models: residual CO2 follows the same law as in the VE models, so 3D
# reference runs and the fine cells of hybrid grids share their fine-scale
# physics.
#
# At each point the residual CO2 scales with the largest saturation reached
# there; only the rest, the free saturation
#   s_free = ((1 − r_w)·s − r_g·s_max)/(1 − r_w − r_g)   (s < s_max; else s),
# flows and holds capillary pressure. The relperm is the drainage curve at
# s_free; the capillary pressure is that of s_free, so pc(0) where only
# residual CO2 is left (co2lab's assumption; JutulDarcy's own hysteresis is
# the alternative, e.g. for sensitivity studies).

"""
    EndpointScalingHysteresis(res_water, res_gas)

co2lab's endpoint-scaling hysteresis as a JutulDarcy hysteresis model for
fine 3D cells: the CO2 relperm is the drainage curve at the free saturation
`((1 − r_w)·s − r_g·s_max)/(1 − r_w − r_g)`. Used as `hysteresis_g` of
JutulDarcy's `ReservoirRelativePermeabilities`.

# Notes

- The relperm needs a drainage and an imbibition CO2 curve. JutulDarcy uses
  the imbibition curve at the largest `s_max = 1 − r_w`, so it must be the
  same law from there: `kr_drain((1 − r_w)·(s − r_g)/(1 − r_w − r_g))`.
"""
struct EndpointScalingHysteresis <: JutulDarcy.AbstractHysteresis
    res_water::Float64  # residual brine saturation
    res_gas::Float64    # residual CO2 saturation
end

JutulDarcy.hysteresis_impl(h::EndpointScalingHysteresis, drain, imb, S, S_max) =
    drain(_free_saturation(S, S_max, h.res_water, h.res_gas))

"""
    FreeSaturationCapillaryPressure(pc, res_water, res_gas)

Capillary pressure of fine 3D cells (CO2 relative to brine, Pa) at the free
CO2 saturation of [`EndpointScalingHysteresis`](@ref): `pc(s_free)`, so
`pc(0)` where only residual CO2 is left, as the VE models assume. Install it
as the model's `CapillaryPressure`.

# Notes

- Needs JutulDarcy's `MaxSaturations` (`add_relperm_parameters!` with
  hysteresis active).
"""
struct FreeSaturationCapillaryPressure{F} <: JutulDarcy.AbstractCapillaryPressure
    pc::F               # capillary pressure (Pa) as a function of CO2 saturation
    res_water::Float64  # residual brine saturation
    res_gas::Float64    # residual CO2 saturation
end

# Modifies: `pc` (one row: CO2 relative to brine).
@jutul_secondary function update_free_saturation_pc!(pc,
                                                    v::FreeSaturationCapillaryPressure,
                                                    model, Saturations, MaxSaturations,
                                                    ix)
    for c in ix
        # the CO2 is the second phase; s_max lags a step (updated before it)
        s = Saturations[2, c]
        pc[1, c] = v.pc(_free_saturation(s, MaxSaturations[2, c], v.res_water,
                                         v.res_gas))
    end
    return pc
end

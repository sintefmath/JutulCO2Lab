# co2lab's endpoint-scaling hysteresis for fine 3D cells (roadmap Stage 6,
# the combined model; design note, point 5): JutulDarcy models whose residual
# CO2 follows the same law as the VE models', so that 3D reference runs share
# their fine-scale physics (and roadmap Stage 12's fine cells can).
#
# At each point the CO2 left behind as residual scales with the largest
# saturation reached there; only the rest, the free saturation
#   s_free = ((1 − r_w)·s − r_g·s_max)/(1 − r_w − r_g)   (s < s_max; else s),
# flows and holds capillary pressure. The relperm is the drainage curve at
# s_free; the capillary pressure is that of s_free, zero where only residual
# CO2 is left (co2lab's assumption: researcher, 2026-10-04, option (c) — this
# for the comparisons, JutulDarcy's own as a sensitivity).

"""
    EndpointScalingHysteresis(res_water, res_gas)

co2lab's endpoint-scaling hysteresis as a JutulDarcy hysteresis model for
fine 3D cells (`hysteresis_g` of JutulDarcy's `ReservoirRelativePermeabilities`,
with a drainage and an imbibition CO2 curve): the CO2 relperm is the drainage
curve at the free saturation `((1 − r_w)·s − r_g·s_max)/(1 − r_w − r_g)`. The
imbibition curve must be the same law from the largest `s_max = 1 − r_w`
(JutulDarcy uses it there): `kr_drain((1 − r_w)·(s − r_g)/(1 − r_w − r_g))`.
Returns a new `EndpointScalingHysteresis`.
"""
struct EndpointScalingHysteresis <: JutulDarcy.AbstractHysteresis
    res_water::Float64  # residual brine saturation
    res_gas::Float64    # residual CO2 saturation
end

JutulDarcy.hysteresis_impl(h::EndpointScalingHysteresis, drain, imb, S, S_max) =
    drain(_free_saturation(S, S_max, h.res_water, h.res_gas))

"""
    FreeSaturationCapillaryPressure(pc, res_water, res_gas)

The capillary pressure of fine 3D cells (CO2 relative to brine, Pa) at the
free CO2 saturation of [`EndpointScalingHysteresis`](@ref): `pc(s_free)`,
`pc` a function of the CO2 saturation. Where only residual CO2 is left it is
`pc(0)`, as the VE models assume. Needs JutulDarcy's `MaxSaturations`
(`add_relperm_parameters!` with hysteresis active). Returns a new
`FreeSaturationCapillaryPressure`; install it as the model's
`CapillaryPressure`.
"""
struct FreeSaturationCapillaryPressure{F} <: JutulDarcy.AbstractCapillaryPressure
    pc::F               # capillary pressure of the CO2 saturation (Pa)
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

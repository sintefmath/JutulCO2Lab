# How the residual CO2 depends on the saturation history (roadmap Stage 12
# step 2, and "Later: co2lab's 'fixed residual' hysteresis"), for VE blocks
# (`LayeredFringe`) and fine 3D cells (`FineCellFluid`) alike.
#
# Each law states, at one depth, the actual CO2 saturation of a column whose
# drainage profile stands at capillary height h (the mobile CO2) after one at
# h_max (the largest reached), and the saturation the CO2 relperm is taken
# at. The column's saturation, relperms and reconstruction integrate these.

"""
    HysteresisLaw

How the residual CO2 depends on the saturation history, for VE blocks and
fine cells alike: [`EndpointScaling`](@ref) (co2lab's default) or
[`FixedResidual`](@ref). Given to [`FineScaleFluid`](@ref) as `hysteresis`,
or to [`LayeredFringe`](@ref).
"""
abstract type HysteresisLaw end

"""
    EndpointScaling()

co2lab's endpoint scaling (its default, `'endpoint scaling'`): at each depth
the residual CO2 scales with the largest saturation reached there,
`r_g·s_max/(1 − r_w)`; the free saturation is `((1 − r_w)·s − r_g·s_max)/(1 −
r_w − r_g)` while `s < s_max`, else `s`. Linear in `(s, s_max)`, so it holds
for column averages too.
"""
struct EndpointScaling <: HysteresisLaw end

"""
    FixedResidual()

A fixed residual saturation, capped by the history (researcher, 2026-10-06):
at each depth the CO2 saturation follows the drainage curve down until it
reaches the residual `min(r_g, s_max)`, and stays there, trapped. In a VE
column at capillary height `h` after `h_max`: the drainage profile of `h`
where it holds at least `r_g`, `r_g` below that in the swept zone, and the
old profile `s_max` where that stayed below `r_g` (all of it trapped). The
CO2 relperm is the drainage curve's where the CO2 is not trapped, else 0.

co2lab's `'fixed residual'` (`addVERelpermCapillaryFringe.m`,
`free_sg_sharp`) puts `r_g` everywhere in the swept zone `[h, h_max]`,
uncapped, even where the drainage profile held less, and none in the fringe
above `h`: more residual CO2 near the old fringe tip than was ever there,
and a height equation that is not monotone for thin plumes (its Newton can
fail; analysis/co2lab-issues.org). The two agree for curves whose fringe
holds at least `r_g` everywhere (a critical saturation of at least `r_g`).
co2lab's reconstruction (`height2finescaleSat.m`) uses this capped rule.

Best with a CO2 relperm curve that is 0 up to `r_g`: one still positive
there drops to 0 when the CO2 is trapped, and Newton needs small steps at
the switch to imbibition (analysis/fixed-residual-hysteresis.md).
"""
struct FixedResidual <: HysteresisLaw end

_free_saturation(::EndpointScaling, s, smax, rw, rg) = _free_saturation(s, smax, rw, rg)

# ----------------------------- Fine cells' rules ------------------------------

# The saturation a fine cell's CO2 relperm (drainage curve) is taken at:
# endpoint scaling's free saturation; with a fixed residual the saturation
# itself, or 0 when trapped (below the maximum and at or below r_g)
_fine_kr_saturation(law::EndpointScaling, s, smax, rw, rg) =
    _free_saturation(law, s, smax, rw, rg)
_fine_kr_saturation(::FixedResidual, s, smax, rw, rg) =
    (Jutul.value(s) < Jutul.value(smax) && Jutul.value(s) <= rg) ? zero(s) : s

# The saturation a fine cell's capillary pressure is taken at: endpoint
# scaling's free saturation; with a fixed residual the saturation itself
# (along the drainage curve; a trapped cell's CO2 does not flow)
_fine_pc_saturation(law::EndpointScaling, s, smax, rw, rg) =
    _free_saturation(law, s, smax, rw, rg)
_fine_pc_saturation(::FixedResidual, s, smax, rw, rg) = s

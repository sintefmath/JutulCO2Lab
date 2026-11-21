# The fluid described at the fine scale: per rock type the CO2 and brine
# relperms and the capillary curve, plus residual saturations, hysteresis law,
# PVT and dissolution. One source for both parts of a hybrid grid, so that VE
# blocks and fine cells agree:
#
# - VE blocks get the upscaled curves: the vertical model `LayeredFringe`, or
#   with no capillarity anywhere its limit `SharpInterfaceIntegrated`;
# - fine 3D cells get the curves themselves (`FineCellFluid`).
#
# co2lab's own VE options (`co2lab_ve_fluid` with keywords) remain for ports
# of co2lab cases. Not part of the fine-scale fluid: co2lab's dissolution rate
# (an upscaled parameter, given to the VE fluid) and the caprock's rugosity
# (the grid's).

# ------------------------------ The fine fluid --------------------------------

"""
    RockCurves(; invPc3D = nothing, kr3D = 3, krw3D = :linear)

One rock type's fine-scale curves, an entry of [`FineScaleFluid`](@ref)'s
`rock_types`.

# Fields

Each is also a constructor keyword of the same name.

- `invPc3D`: brine saturation as a function of capillary pressure (Pa);
  `nothing`: no capillarity.
- `kr3D`: CO2 relperm as a function of CO2 saturation, or a Corey exponent
  `β` for `s^β`.
- `krw3D`: brine relperm as a function of brine saturation, or `:linear` for
  `(s_w − r_w)/(1 − r_w)`.

# Notes

- **Performance:** the curves are evaluated in the simulator's innermost
  loops (the capillary fringe's quadrature, every Newton iteration, with AD
  numbers). A curve that reads a non-constant global variable, e.g. a
  parameter from a script's `settings`, is type-unstable and can make a run
  many times slower. Give the parameters to the curve as typed values
  instead: build it inside a function, or in a `let` block, so that the
  closure captures them, e.g.

      kr3D = let n = settings.n_g
          s -> s^n
      end
"""
struct RockCurves{P, K, W}
    invPc3D::P  # brine saturation of capillary pressure (Pa), or nothing
    kr3D::K     # CO2 relperm of CO2 saturation, or a Corey exponent
    krw3D::W    # brine relperm of brine saturation, or :linear
end
RockCurves(; invPc3D = nothing, kr3D = 3, krw3D = :linear) =
    RockCurves(invPc3D, kr3D, krw3D)

"""
    FineScaleFluid(; rock_types = [RockCurves()], residual = (0.0, 0.0),
                     hysteresis = EndpointScaling(), fine_hysteresis = :law,
                     kwarg...)

The fluid described at the fine scale, the single source for a hybrid grid's
VE blocks and fine cells. Turn it into a VE fluid with
[`co2lab_ve_fluid`](@ref)`(fluid, mesh, columns; …)`.

# Fields

- `rock_types`: the curves of each rock type ([`RockCurves`](@ref)).
- `residual`: `(r_w, r_g)`, residual brine and CO2 saturations, for all rock
  types.
- `hysteresis`: the hysteresis law ([`HysteresisLaw`](@ref)) of VE blocks and
  fine cells.
- `fine_hysteresis`: the fine cells' hysteresis, `:law` or `:killough` (see
  Keywords).
- `options`: the PVT and dissolution options of [`co2lab_ve_fluid`](@ref),
  as given in `kwarg...`.

# Keywords

- `rock_types`: the curves of each rock type ([`RockCurves`](@ref)).
- `residual`: `(r_w, r_g)`, residual brine and CO2 saturations, for all rock
  types.
- `hysteresis`: the hysteresis law ([`HysteresisLaw`](@ref)) of VE blocks and
  fine cells.
- `fine_hysteresis`: the fine cells' hysteresis: `:law` (the law above), or
  `:killough` (JutulDarcy's Killough model, with the capillary pressure of the
  total saturation; VE blocks keep the law, and a warning notes the mismatch).
- `kwarg...`: the PVT and dissolution options of [`co2lab_ve_fluid`](@ref)
  (`co2_rho_pvt`, …, `dissolution`, `dis_max`, `dis_max_mass`). Throws for
  `vertical`, `krmax` and `dis_rate`, which are not fine-scale options.
"""
struct FineScaleFluid{L <: HysteresisLaw}
    rock_types::Vector{RockCurves}  # fine-scale curves per rock type
    residual::NTuple{2, Float64}    # residual brine and CO2 saturations
    hysteresis::L                   # hysteresis law (VE blocks and fine cells)
    fine_hysteresis::Symbol         # fine cells: :law, or :killough (JutulDarcy's)
    options::NamedTuple             # co2lab_ve_fluid's PVT and dissolution options
end

function FineScaleFluid(; rock_types = [RockCurves()], residual = (0.0, 0.0),
                        hysteresis::HysteresisLaw = EndpointScaling(),
                        fine_hysteresis::Symbol = :law, kwarg...)
    fine_hysteresis in (:law, :killough) ||
        throw(ArgumentError("fine_hysteresis: :law or :killough"))
    bad = intersect(keys(kwarg), (:vertical, :krmax, :dis_rate))
    isempty(bad) || throw(ArgumentError("not a fine-scale option: $(first(bad))"))
    return FineScaleFluid(collect(RockCurves, rock_types), Float64.(Tuple(residual)),
                          hysteresis, fine_hysteresis, NamedTuple(kwarg))
end

"""
    co2lab_ve_fluid(fluid::FineScaleFluid, mesh, columns; poro, perm_x,
                    perm_y = nothing, rock_type = nothing, dis_rate = 0.0)
        -> CO2LabVEFluid

The VE fluid of a fine-scale `fluid` on the blocks of `columns`: its curves
upscaled over each column, with residual, PVT and dissolution taken from
`fluid`.

# Arguments

- `fluid`: the fine-scale fluid ([`FineScaleFluid`](@ref)).
- `mesh`, `columns`: the 3D mesh and its columns ([`ve_columns`](@ref)).
- `poro`, `perm_x`: porosity and x-permeability (m²) per cell of `mesh`, as
  given to [`ve_domain`](@ref) (net rock with net-to-gross: [`net_rock`](@ref)).
- `perm_y`: y-permeability per cell, only checked against `perm_x`.
- `rock_type`: rock type per cell of `mesh` (index into
  `fluid.rock_types`); `nothing`: rock type 1 everywhere.
- `dis_rate`: co2lab's upscaled dissolution rate, for VE blocks (fine cells
  dissolve instantaneously).

# Returns

A new `CO2LabVEFluid` whose vertical model is
- [`LayeredFringe`](@ref) of the curves and hysteresis law, or,
- with no capillarity in any rock type, [`SharpInterfaceIntegrated`](@ref)
  (the same limit, where both laws agree) with the curves' endpoints; the
  rock types must then share their relperm endpoints.

Its field `fine` carries `(fluid, rock_type)` for the hybrid grid's fine
cells.
"""
function co2lab_ve_fluid(fluid::FineScaleFluid, mesh, columns::VEColumns; poro, perm_x,
                         perm_y = nothing, rock_type = nothing, dis_rate = 0.0)
    types = collect(Int, something(rock_type, ones(Int, number_of_cells(mesh))))
    maximum(types) <= length(fluid.rock_types) ||
        throw(ArgumentError("rock_type beyond the fluid's $(length(fluid.rock_types))"))
    rc, (rw, rg) = fluid.rock_types, fluid.residual
    fluid.fine_hysteresis == :killough && rg > 0 &&
        @warn "Killough's hysteresis in fine cells, $(fluid.hysteresis) in VE blocks"
    if all(r -> isnothing(r.invPc3D), rc)
        vertical = SharpInterfaceIntegrated(mesh, columns; poro, perm_x, perm_y)
        krmax = _sharp_endpoints(rc, rw, rg)
    else
        vertical = LayeredFringe(mesh, columns; poro, perm_x, perm_y, rock_type = types,
                                 invPc3D = [r.invPc3D for r in rc],
                                 kr3D = [r.kr3D for r in rc],
                                 krw3D = [r.krw3D for r in rc],
                                 hysteresis = fluid.hysteresis)
        krmax = (-1.0, -1.0)
    end
    ve = co2lab_ve_fluid(; vertical, residual = (rw, rg), krmax, dis_rate,
                         fluid.options...)
    fields = (getfield(ve, k) for k in fieldnames(CO2LabVEFluid)[1:end-1])
    return CO2LabVEFluid(fields..., (fluid = fluid, rock_type = types))
end

# The sharp interface's relperm endpoints from the curves: CO2 at 1 − r_w,
# brine (residual zone) at 1 − r_g; throws unless all rock types agree
function _sharp_endpoints(rc, rw, rg)
    ends = unique((_brine_curve(r.krw3D, rw)(1 - rg), _co2_curve(r.kr3D)(1 - rw))
                  for r in rc)
    length(ends) == 1 ||
        throw(ArgumentError("without capillarity the rock types need the same " *
                            "relperm endpoints (SharpInterfaceIntegrated has one pair)"))
    return only(ends)
end

"""
    _saturation_kinks(fluid::FineScaleFluid) -> Vector{Vector{Float64}}

Per rock type, the CO2 saturations where a fine cell's functions kink, at
which [`FineGasSaturation`](@ref) stops Newton updates: the critical
saturation (the largest with zero CO2 relperm, if above 0) and the residual
CO2 `r_g` (if above 0, where the hysteresis laws switch). Returns new sorted
vectors, coinciding kinks merged.
"""
function _saturation_kinks(fluid::FineScaleFluid)
    rw, rg = fluid.residual
    return map(fluid.rock_types) do r
        kr = _co2_curve(r.kr3D)
        ks = Float64[]
        # the critical saturation: kr zero below, positive above
        if kr(1e-12) == 0 && kr(1 - rw) > 0
            lo, hi = 0.0, 1 - rw
            for _ in 1:60
                mid = (lo + hi) / 2
                kr(mid) > 0 ? (hi = mid) : (lo = mid)
            end
            lo > 0 && push!(ks, lo)
        end
        rg > 0 && push!(ks, rg)
        # one kink where the two coincide (to round-off)
        sort!(ks)
        ks[[i == 1 || ks[i] - ks[i-1] > 1e-12 for i in eachindex(ks)]]
    end
end

# The curves as functions
_co2_curve(kr3D::Number) = s -> max(s, 0.0)^kr3D
_co2_curve(kr3D) = kr3D
_brine_curve(::Symbol, rw) = sw -> max(sw - rw, 0.0) / (1 - rw)
_brine_curve(krw3D, rw) = krw3D

# -------------------------------- Fine cells ----------------------------------

"""
    FineCellFluid

The constitutive functions of a hybrid grid's fine 3D cells, built from a
[`FineScaleFluid`](@ref) and each cell's rock type: per rock type, the
relperms and capillary pressure as tables in the saturation. The imbibition
curve (from the largest maximum, 1 − r_w) is for Killough's model.

# Fields

- `drain`: per rock type, the CO2 relperm table of the drainage curve.
- `imb`: per rock type, the CO2 relperm table of the imbibition curve from
  the largest maximum CO2 saturation, 1 − r_w.
- `brine`: per rock type, the brine relperm table in the brine saturation.
- `pc`: per rock type, the capillary pressure table in the CO2 saturation
  (Pa).
- `rock_type`: per domain cell, its rock type; 0 for a VE block.
- `law`: the hysteresis law ([`HysteresisLaw`](@ref)).
- `killough`: `true` for JutulDarcy's Killough hysteresis instead of `law`.
- `res_water`: residual brine saturation.
- `res_gas`: residual CO2 saturation.
"""
struct FineCellFluid{K, I, L}
    drain::Vector{K}        # per rock type: CO2 relperm, drainage
    imb::Vector{K}          # per rock type: CO2 relperm, imbibition from 1 − r_w
    brine::Vector{I}        # per rock type: brine relperm of brine saturation
    pc::Vector{I}           # per rock type: capillary pressure of CO2 saturation (Pa)
    rock_type::Vector{Int}  # per domain cell (0: VE block)
    law::L                  # hysteresis law
    killough::Bool          # Killough's hysteresis (JutulDarcy's) instead of the law
    res_water::Float64      # residual brine saturation
    res_gas::Float64        # residual CO2 saturation
end

"""
    _fine_cell_fluid(fine, domain) -> FineCellFluid

The fine cells' functions for `domain` (a hybrid `ve_domain`) from a VE
fluid's `fine` (`(fluid, rock_type)`, rock types per cell of the 3D mesh).
Returns a new `FineCellFluid`.
"""
function _fine_cell_fluid(fine, domain)
    ff, types = fine.fluid, fine.rock_type
    (rw, rg), rc = ff.residual, ff.rock_types
    coarse = Jutul.physical_representation(domain)
    is_ve = domain[:is_ve]
    rt = [is_ve[b] ? 0 : types[coarse.partition_lookup[b][1]] for b in eachindex(is_ve)]
    s = collect(range(0.0, 1.0; length = 2001))
    tab(f) = JutulDarcy.PhaseRelativePermeability(s, f.(s); label = :g)
    interp(f) = Jutul.get_1d_interpolator(s, f.(s))
    # CO2: drainage up to 1 − r_w; imbibition from there (endpoint scaling)
    co2(r) = x -> _co2_curve(r.kr3D)(min(x, 1 - rw))
    imb(r) = rg > 0 ? (x -> co2(r)((1 - rw) * max(x - rg, 0.0) / (1 - rw - rg))) : co2(r)
    return FineCellFluid([tab(co2(r)) for r in rc], [tab(imb(r)) for r in rc],
                         [interp(_brine_curve(r.krw3D, rw)) for r in rc],
                         [interp(_capillary_of_co2(r.invPc3D, rw)) for r in rc], rt,
                         ff.hysteresis, ff.fine_hysteresis == :killough && rg > 0,
                         rw, rg)
end

# The capillary pressure as a function of CO2 saturation (the inverse of
# invPc3D); at zero CO2 the entry pressure, the limit from above
_capillary_of_co2(::Nothing, rw) = sg -> 0.0
function _capillary_of_co2(invpc, rw)
    sgp(p) = 1 - invpc(max(p, 0.0))
    pmax = _upper_capillary_pressure(invpc, rw, 0.0)
    pe = _entry_pressure(sgp, pmax)
    return function (sg)
        sg <= 0 && return pe
        sg >= sgp(pmax) && return pmax
        lo, hi = pe, pmax
        for _ in 1:100
            mid = (lo + hi) / 2
            sgp(mid) < sg ? (lo = mid) : (hi = mid)
        end
        return hi
    end
end

"""
    _fine_relperm(f, i, sG, smax) -> (krW, krG)

Fine cell `i`'s relperms (brine, CO2) at CO2 saturation `sG` and maximum
`smax`. CO2: the drainage curve at the hysteresis law's saturation (the
free one, or `sG` itself unless trapped), or Killough's model. Brine: its
curve at `1 − sG`.
"""
function _fine_relperm(f::FineCellFluid, i, sG, smax)
    t = f.rock_type[i]
    krg = f.killough ?
          JutulDarcy.kr_hysteresis(JutulDarcy.KilloughHysteresis(), f.drain[t], f.imb[t],
                                   sG, smax) :
          f.drain[t](_fine_kr_saturation(f.law, sG, smax, f.res_water, f.res_gas))
    return (f.brine[t](1 - sG), krg)
end

"""
    _fine_capillary_pressure(f, i, sG, smax) -> Real

Fine cell `i`'s capillary pressure (CO2 relative to brine), evaluated at
- the hysteresis law's saturation: the free one for endpoint scaling (zero
  free CO2: the entry pressure), `sG` for a fixed residual;
- `sG` with Killough's model.
"""
function _fine_capillary_pressure(f::FineCellFluid, i, sG, smax)
    t = f.rock_type[i]
    s = f.killough ? sG : _fine_pc_saturation(f.law, sG, smax, f.res_water, f.res_gas)
    return f.pc[t](s)
end

# The capillary fringe in a column whose rock varies with depth (roadmap Stage
# 6, the combined model; design: analysis/stage6-upscaling-proposal.md). Not
# in co2lab, which has the fringe ('P-scaled table') and the layers
# ('sharp_interface_integrated') only one at a time.
#
# The fine-scale capillary pressure is linear in depth, pc(z) = Δρ g (h − z),
# h the capillary height. Within one fine cell the rock is constant, so the
# column's integrals split into one P-type term per cell, from the fringe
# table's integrals over capillary pressure (SP, krP; `PScaledTable`):
#   ∫ s dz over [a, b] = (SP(Δρ g (h − a)) − SP(Δρ g (h − b))) / (Δρ g),
# weighted by the cell's pore volume per depth (saturation) or permeability
# (relperm). The current Δρ enters every evaluation; nothing is stored per
# column but the rock. h is solved from the free saturation on plain values,
# its derivatives taken by one Newton step with the AD values.
#
# Decided with the researcher (2026-10-04; the design note's points 1–7):
# - hysteresis by co2lab's endpoint scaling, linear in s and s_max at each
#   depth, so it holds for the column averages too: h from the free
#   saturation, h_max from sGmax, both at the current Δρ (an approximation:
#   the maximum profile was set at the Δρ of its time);
# - the brine relperm integrated from a fine-scale curve `krw3D` over the
#   actual saturation (free and residual), linear by default — co2lab's
#   krW_simple in drainage, less by r_w/(1 − r_w)·(S − S_free) with residual CO2;
# - physical capillary curves only (`co2lab_capillary_curve` converts
#   co2lab's dimensionless (C, alpha));
# - no caprock rugosity (a second pass);
# - built on integrals over any depth interval, which roadmap Stage 12's
#   transition flux needs.
#
# Uniform rock gives `CapillaryFringeTable`'s heights and CO2 relperm.

# --------------------------------- Settings -----------------------------------

"""
    LayeredFringe(mesh, columns; poro, perm_x, invPc3D, kr3D = 3, krw3D = :linear,
                  perm_y = nothing, geometry = ve_column_geometry(mesh, columns),
                  tolerance = 1e-6, tip_saturation = 1e-6,
                  touched_saturation = 1e-3) -> LayeredFringe

A capillary fringe in VE columns whose rock varies with depth: the columns
`columns` of the 3D `mesh`, with the fine cells' porosity `poro` and
x-permeability `perm_x` (per cell of `mesh`, as given to
[`ve_domain`](@ref); `perm_y` only checked, as for
[`SharpInterfaceIntegrated`](@ref)). One fine-scale capillary curve for all
cells:
- `invPc3D`: the fine-scale brine saturation as a function of capillary
  pressure (Pa); a physical curve (for co2lab's dimensionless `(C, alpha)`
  see [`co2lab_capillary_curve`](@ref));
- `kr3D`: the fine-scale CO2 relperm of the CO2 saturation, or an exponent
  `β` for `s^β`;
- `krw3D`: the fine-scale brine relperm of the brine saturation, or
  `:linear` for `(s_w − r_w)/(1 − r_w)`;
- `tolerance`, `tip_saturation`, `touched_saturation`: as for
  [`CapillaryFringeTable`](@ref).

Saturations are integrated by pore volume, relperms by permeability, over
the fringe's profile cell by cell; hysteresis is co2lab's endpoint scaling.
Does not combine with the caprock's rugosity: [`setup_co2lab_model`](@ref)
refuses a domain with a non-zero `:rugosity`.

Returns a new `LayeredFringe` to give to [`co2lab_ve_fluid`](@ref) as
`vertical`; it belongs to this grid. The arguments are only read.
"""
struct LayeredFringe{F, K, W} <: VerticalModel
    cells::SharpInterfaceIntegrated     # per column: cell boundaries, ∫φ·area, ∫k_x dz
    fringe::CapillaryFringeTable{F, K}  # fine-scale curves and table settings
    krw3D::W                            # brine relperm of brine saturation, or :linear
end

function LayeredFringe(mesh, columns::VEColumns; poro, perm_x, invPc3D, kr3D = 3,
                       krw3D = :linear, perm_y = nothing,
                       geometry = ve_column_geometry(mesh, columns),
                       tolerance = 1e-6, tip_saturation = 1e-6,
                       touched_saturation = 1e-3)
    invPc3D isa Tuple &&
        throw(ArgumentError("LayeredFringe takes a physical capillary curve; " *
                            "co2lab_capillary_curve converts co2lab's (C, alpha)"))
    krw3D isa Symbol && krw3D != :linear &&
        throw(ArgumentError("krw3D: a function of the brine saturation, or :linear"))
    cells = SharpInterfaceIntegrated(mesh, columns; poro, perm_x, perm_y, geometry)
    fringe = CapillaryFringeTable(; invPc3D, kr3D, tolerance, tip_saturation,
                                  touched_saturation)
    return LayeredFringe(cells, fringe, krw3D)
end

"""
    co2lab_capillary_curve(C, alpha; height, drho, res_water = 0.0) -> Function

co2lab's dimensionless fine-scale capillary curve `(C, alpha)` as a physical
one: the brine saturation `s_w(p) = max((C'/(p + C'))^(1/alpha), res_water)`
of the capillary pressure `p` (Pa), with `C' = C·g·height·drho`. co2lab takes
`height` as the grid's tallest column and `drho` as its reference density
difference ([`CapillaryFringeTable`](@ref) does so); here they are given.
Returns a new function.
"""
function co2lab_capillary_curve(C, alpha; height, drho, res_water = 0.0)
    Cs = C * Jutul.gravity_constant * height * drho
    return p -> max((Cs / (max(p, 0.0) + Cs))^(1 / alpha), res_water)
end

"""
    LayeredFringeModel

The built form of [`LayeredFringe`](@ref) for a domain and fluid: the
columns' rock, the fringe table and the residual saturations.
"""
struct LayeredFringeModel{T, W} <: VerticalModel
    cells::SharpInterfaceIntegrated  # per column: cell boundaries, ∫φ·area, ∫k_x dz
    table::T                         # fine-scale integrals over capillary pressure
    krw3D::W                         # brine relperm of brine saturation, or :linear
    res_water::Float64               # residual brine saturation
    res_gas::Float64                 # residual CO2 saturation
end

# The model belongs to its grid; the table as CapillaryFringeTable's
function _vertical_model(spec::LayeredFringe, domain, fluid)
    _refuse_rugosity(domain, "LayeredFringe")
    _vertical_model(spec.cells, domain, fluid)
    table = _vertical_model(spec.fringe, domain, fluid)
    return LayeredFringeModel(spec.cells, table, spec.krw3D, fluid.res_water,
                              fluid.res_gas)
end

# ----------------------------- Interval integrals -----------------------------

# Column i's cells, as indices a into m.cells (the cell from z[a] to z[a+1])
_cell_range(c::SharpInterfaceIntegrated, i) = c.pos[i]:(c.pos[i+1] - 2)

# ∫ of table Y (SP or krP) over depths [a, b] for capillary height h: the
# integral of the fine-scale quantity over that depth range (dg = Δρ g)
function _depth_integral(t, Y, dY, h, a, b, dg)
    return (_table_value(t, Y, dY, dg * (h - a)) - _table_value(t, Y, dY, dg * (h - b))) /
           dg
end

"""
    _interval_integrals(m, i, T, B, h, h_max, drho) -> NamedTuple

Column `i`'s integrals over the depths `[T, B]` below its caprock, for the
capillary heights `h` (free CO2) and `h_max` (maximum profile) and density
difference `drho`, cell by cell: `free`, `max`, `co2` the CO2 saturation of
the free, maximum and actual profiles times the pore volume (m³; actual =
`(1 − C)·free + C·max`, `C = r_g/(1 − r_w)`, endpoint scaling), `kfree`,
`kmax` the free and maximum profiles' CO2 saturation times `k_x dz`, `krG`,
`krW` the CO2 and brine relperms times `k_x dz` (m³), `pv` and `kz` the
interval's pore volume and `∫ k_x dz`. The column's quantities are those over
`[0, H]`; roadmap Stage 12's transition flux asks for others. Arguments may
carry AD derivatives. Returns a new NamedTuple.
"""
function _interval_integrals(m::LayeredFringeModel, i, T, B, h, h_max, drho)
    c, t = m.cells, m.table
    dg = drho * Jutul.gravity_constant
    C = m.res_gas / (1 - m.res_water)
    z0 = zero(h) + zero(h_max) + zero(drho)
    free, smax, co2, kf, km, krG, krW = z0, z0, z0, z0, z0, z0, z0
    pv, kz = 0.0, 0.0
    hv = max(Jutul.value(h), Jutul.value(h_max))
    for a in _cell_range(c, i)
        top, bot = max(c.z[a], T), min(c.z[a+1], B)
        bot > top || continue
        dz = c.z[a+1] - c.z[a]
        phiA, k = (c.V[a+1] - c.V[a]) / dz, (c.K[a+1] - c.K[a]) / dz
        pv, kz = pv + phiA * (bot - top), kz + k * (bot - top)
        # below both profiles: brine only (the top cell of an empty column
        # counts, for the slope at S = 0)
        if top > hv
            krW += k * (bot - top) * _brine_relperm_full(m)
            continue
        end

        sf = _depth_integral(t, t.SP, t.s, h, top, bot, dg)
        sm = _depth_integral(t, t.SP, t.s, h_max, top, bot, dg)
        s = (1 - C) * sf + C * sm
        free, smax, co2 = free + phiA * sf, smax + phiA * sm, co2 + phiA * s
        kf, km = kf + k * sf, km + k * sm
        krG += k * _depth_integral(t, t.krP, t.kr, h, top, bot, dg)
        krW += k * _brine_integral(m, top, bot, s, h, h_max, dg)
    end
    return (free = free, max = smax, co2 = co2, kfree = kf, kmax = km, krG = krG,
            krW = krW, pv = pv, kz = kz)
end

# The brine relperm of brine alone (s_w = 1)
_brine_relperm_full(m::LayeredFringeModel{<:Any, Symbol}) = 1.0
_brine_relperm_full(m::LayeredFringeModel) = m.krw3D(1.0)

# ∫ krw3D(1 − s(z)) dz over [top, bot] (`s` the actual profile's ∫ s dz there):
# linear, (s_w − r_w)/(1 − r_w), in closed form
_brine_integral(m::LayeredFringeModel{<:Any, Symbol}, top, bot, s, h, h_max, dg) =
    (bot - top) - s / (1 - m.res_water)

# A function krw3D: Gauss–Legendre between the profiles' kinks (where the
# capillary pressure of either reaches 0 or the entry pressure). The kinks
# carry their AD derivatives: the saturation may jump there (an entry
# pressure, no capillarity), and the moving jump is part of the derivative.
function _brine_integral(m::LayeredFringeModel, top, bot, s, h, h_max, dg)
    t = m.table
    C = m.res_gas / (1 - m.res_water)
    sg(p) = 1 - t.invpc(_max_c(p, 0.0))
    sw(z) = 1 - (1 - C) * sg(dg * (h - z)) - C * sg(dg * (h_max - z))
    z0 = zero(h) + zero(h_max) + zero(dg)
    kinks = [z0 + x - e / dg for x in (h, h_max) for e in (0.0, t.entry)]
    inside = filter(z -> top < Jutul.value(z) < bot, kinks)
    cuts = sort!([z0 + top; inside; z0 + bot]; by = Jutul.value)
    return sum(_gl_integral(z -> m.krw3D(sw(z)), a, b)
               for (a, b) in zip(cuts[1:end-1], cuts[2:end]))
end

# ---------------------------------- Heights -----------------------------------

# Column i's height and total pore volume (its cells')
_column_extent(c::SharpInterfaceIntegrated, i) = (c.z[c.pos[i+1] - 1], _total(c, i, c.V))

"""
    _layered_height(m, i, S, drho) -> Real

The capillary height of column `i` whose drainage profile holds the CO2
saturation `S` (pore-volume weighted), at density difference `drho`: the root
of `Σ φA ∫ s dz − S·V` over the column. Solved on plain values to round-off,
the derivatives by one Newton step with the AD values; below the tip
saturation linear in `S`, as [`CapillaryFringeTable`](@ref). Zero for
`drho ≤ 0`.
"""
function _layered_height(m::LayeredFringeModel, i, S, drho)
    # no buoyancy (CO2 as dense as brine, an unphysical state Newton may pass
    # through): no height, as CapillaryFringeTable's
    Jutul.value(drho) > 0 || return zero(S) + zero(drho)
    tip = m.table.tip
    Jutul.value(S) < tip && return _layered_newton(m, i, tip, drho) * S / tip
    return _layered_newton(m, i, S, drho)
end

# The root on plain values, then one Newton step with the AD values
function _layered_newton(m, i, S, drho)
    _, V = _column_extent(m.cells, i)
    h = _solve_layered(m, i, Jutul.value(S) * V, Jutul.value(drho))
    F = _profile_volume(m, i, h, drho) - S * V
    dF = _profile_slope(m, i, h, Jutul.value(drho))
    # a full column: h is not determined (any larger one fills it too)
    dF > 0 || return zero(F) + h
    return h - F / dF
end

# The drainage profile's pore volume of CO2 at height h, Σ φA ∫ s dz (the
# cells below h hold none); h and drho may carry AD derivatives
function _profile_volume(m, i, h, drho)
    c, t = m.cells, m.table
    dg = drho * Jutul.gravity_constant
    v = zero(h) + zero(dg)
    for a in _cell_range(c, i)
        c.z[a] > Jutul.value(h) && break
        dz = c.z[a+1] - c.z[a]
        dz > 0 || continue
        phiA = (c.V[a+1] - c.V[a]) / dz
        v += phiA * _depth_integral(t, t.SP, t.s, h, c.z[a], c.z[a+1], dg)
    end
    return v
end

# d/dh of the drainage profile's pore volume of CO2: Σ φA (s(top) − s(bottom))
function _profile_slope(m, i, h, drho)
    c, t = m.cells, m.table
    dg = drho * Jutul.gravity_constant
    slope(x) = _table_slope(t, t.SP, t.s, x)
    d = 0.0
    for a in _cell_range(c, i)
        c.z[a] > h && break
        dz = c.z[a+1] - c.z[a]
        dz > 0 || continue
        phiA = (c.V[a+1] - c.V[a]) / dz
        d += phiA * (slope(dg * (h - c.z[a])) - slope(dg * (h - c.z[a+1])))
    end
    return d
end

# h with the drainage profile holding the pore volume `v` of CO2 (plain
# values); Newton safeguarded by bisection, to round-off. Beyond the table's
# end the whole column is at its saturation: more CO2 gives that height.
function _solve_layered(m, i, v, drho)
    v <= 0 && return 0.0
    H, V = _column_extent(m.cells, i)
    f(h) = _profile_volume(m, i, h, drho) - v
    lo, hi = 0.0, m.table.p[end] / (drho * Jutul.gravity_constant) + H
    f(hi) <= 0 && return hi
    # start near the root: the plume fully saturated (1 − r_w) at the
    # column's mean porosity
    h = clamp(v / (V / H * (1 - m.res_water)), lo, hi)
    for _ in 1:200
        r = f(h)
        r > 0 ? (hi = h) : (lo = h)
        abs(r) <= 1e-13 * max(v, eps()) && return h
        d = _profile_slope(m, i, h, drho)
        next = d > 0 ? h - r / d : (lo + hi) / 2
        h = lo < next < hi ? next : (lo + hi) / 2
        hi - lo <= 4eps(hi) && return h
    end
    return h
end

# ------------------------- Vertical-model interface ---------------------------

# Capillary heights of the free saturation and of the maximum saturation
function _column_heights(m::LayeredFringeModel, i, sG, sGmax, H, drho, rw, rg)
    h = _layered_height(m, i, _free_saturation(sG, sGmax, rw, rg), drho)
    return (h, _layered_height(m, i, sGmax, drho))
end

# Both relperms integrated over the column, weighted by permeability; the
# endpoints krw, krg are the fine-scale curves' (unused)
function _column_relperm(m::LayeredFringeModel, i, sG, sGmax, h, h_max, H, drho, dh,
                         krw, krg, rw, rg)
    Hc, _ = _column_extent(m.cells, i)
    I = _interval_integrals(m, i, 0.0, Hc, h, h_max, drho)
    Sf = _free_saturation(sG, sGmax, rw, rg)
    krW = _tip_brine_relperm(m, i, Sf, sGmax, h_max, drho, Hc, I)
    return (krW / I.kz, I.krG / I.kz)
end

# The brine relperm (times ∫ k dz) below the tip saturation. There h is linear
# in S, so the profile at h holds CO2 ∝ S², not S, and the integrated brine
# relperm would have no slope at S = 0. Its deficit (∫ k dz less it) is taken
# linear in each profile's saturation instead, exact at the tip: for the
# linear krw3D the deficit is linear in the profiles' CO2, so this is exact
# (co2lab's krW_simple in a homogeneous column); for another krw3D to first
# order.
function _tip_brine_relperm(m::LayeredFringeModel, i, Sf, Sm, h_max, drho, Hc, I)
    tip = m.table.tip
    sf, sm = Jutul.value(Sf), Jutul.value(Sm)
    (sf >= tip && sm >= tip) && return I.krW
    ht = _layered_newton(m, i, tip, drho)
    at(hf, hm) = _interval_integrals(m, i, 0.0, Hc, hf, hm, drho)
    if m.krw3D isa Symbol
        # linear krw3D: the deficit is ((1 − C)·kfree + C·kmax)/(1 − r_w)
        C = m.res_gas / (1 - m.res_water)
        Kt = at(ht, ht).kfree
        Kf = sf < tip ? Kt * Sf / tip : I.kfree
        Km = sm < tip ? Kt * Sm / tip : I.kmax
        return I.kz - ((1 - C) * Kf + C * Km) / (1 - m.res_water)
    elseif sm < tip
        # the whole column below the tip: along the ray to S_max = tip (an
        # empty column along the drainage ray, keeping the slope at S = 0);
        # for a numerically empty column the ratio as a plain value, whose
        # derivatives would overflow (they are scaled by S_max/tip anyway)
        ratio = sm > 1e-12 * tip ? Sf / Sm : (sm > 0 ? sf / sm : 1.0)
        # (below the tip h is linear in S, as everywhere: no solve at tiny S)
        J = at(_layered_height(m, i, tip * ratio, drho), ht)
        return I.kz - (I.kz - J.krW) * Sm / tip
    end
    # only the free profile below the tip
    D0, Dt = I.kz - at(0.0, h_max).krW, I.kz - at(ht, h_max).krW
    return I.kz - (D0 + (Dt - D0) * Sf / tip)
end

# Each fine cell's average of the actual and the maximum profile (cell
# bottoms `zb`; the rock is constant within a cell, so plain depth averages)
function _fine_saturations(m::LayeredFringeModel, i, zb, h, h_max, drho, rw, rg)
    t = m.table
    dg = drho * Jutul.gravity_constant
    C = rg / (1 - rw)
    tops = [0.0; zb[1:end-1]]
    avg(hh, a, b) = b > a ? _depth_integral(t, t.SP, t.s, hh, a, b, dg) / (b - a) : 0.0
    s_max = [avg(h_max, a, b) for (a, b) in zip(tops, zb)]
    s_free = [avg(h, a, b) for (a, b) in zip(tops, zb)]
    return ((1 - C) .* s_free .+ C .* s_max, s_max)
end

# The pore volume below the touched depth (h_max less the entry height; below
# the touched saturation linear in smax, as CapillaryFringeTable's)
function _untouched_brine(m::LayeredFringeModel, i, smax, h_max, H, drho)
    t, c = m.table, m.cells
    dg = drho * Jutul.gravity_constant
    if Jutul.value(smax) < t.touched
        h_t = _max_c(_layered_height(m, i, t.touched, drho) - t.entry / dg, 0.0) *
              smax / t.touched
    else
        h_t = _max_c(h_max - t.entry / dg, 0.0)
    end
    Hc, V = _column_extent(c, i)
    return (V - _integral(c, i, c.V, _min_c(h_t, Hc))) / V
end

# The inventory's split as CapillaryFringeTable's, each profile's fraction
# between two depths by pore volume; co2lab_compatible: co2lab's by height
function _category_saturations(m::LayeredFringeModel, S, Smax, h, h_max, col, drho,
                               rw, rg; co2lab_compatible::Bool = false)
    co2lab_compatible && return _height_split(S, Smax, h, h_max, col, rw, rg)
    S_immob = Smax .* rg ./ (1 - rw)
    S_mob = S .- S_immob
    # the fraction of a drainage profile (height hh) between depths a and b
    function frac(i, hh, a, b)
        Hc, _ = _column_extent(m.cells, i)
        total = _interval_integrals(m, i, 0.0, Hc, hh, hh, drho[i]).free
        total > 0 || return 0.0
        part = _interval_integrals(m, i, clamp(a, 0.0, Hc), clamp(b, 0.0, Hc), hh, hh,
                                   drho[i]).free
        return max(part, 0.0) / total
    end
    mob(i, a, b) = frac(i, h[i], a, b)
    imm(i, a, b) = frac(i, h_max[i], a, b)
    n, dh, zt = eachindex(S), col.dh, col.zt
    struct_end(i) = max(zt[i], dh[i])
    immob_sub = [S_immob[i] * imm(i, 0.0, dh[i]) for i in n]
    mobile_sub = [S_mob[i] * mob(i, 0.0, dh[i]) for i in n]
    immob_struct = [S_immob[i] * imm(i, dh[i], struct_end(i)) for i in n]
    mobile_struct = [S_mob[i] * mob(i, dh[i], struct_end(i)) for i in n]
    # residual below the plume: the immobile part below h and the trap
    res = [S_immob[i] * imm(i, max(h[i], zt[i]), col.H[i]) for i in n]
    return (immob_sub = immob_sub, mobile_sub = mobile_sub,
            immob_struct = immob_struct, mobile_struct = mobile_struct,
            immob_free = S_immob .- immob_struct .- immob_sub,
            mobile_free = S_mob .- mobile_struct .- mobile_sub, res = res)
end

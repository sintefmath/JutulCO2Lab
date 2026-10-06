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
#   the maximum profile was set at the Δρ of its time); or by a fixed
#   residual, capped by the history (`FixedResidual`; researcher,
#   2026-10-06): h then solves the column's CO2 (monotone in h), with the
#   depths where either profile holds r_g as breakpoints in each cell;
# - the brine relperm integrated from a fine-scale curve `krw3D` over the
#   actual saturation (free and residual), linear by default — co2lab's
#   krW_simple in drainage, less by r_w/(1 − r_w)·(S − S_free) with residual CO2;
# - physical capillary curves only (`co2lab_capillary_curve` converts
#   co2lab's dimensionless (C, alpha));
# - the caprock's rugosity dh: the CO2 filling the subscale relief (the top dh
#   of the column, at the fringe's saturation there) is immobile, and the
#   relief takes no part in lateral flow, brine included — both relperms are
#   integrated from dh (co2lab's accretion layer of the sharp interface;
#   researcher, 2026-10-04: the second pass);
# - built on integrals over any depth interval, which roadmap Stage 12's
#   transition flux needs.
#
# Uniform rock gives `CapillaryFringeTable`'s heights and CO2 relperm.

# --------------------------------- Settings -----------------------------------

"""
    LayeredFringe(mesh, columns; poro, perm_x, invPc3D, kr3D = 3, krw3D = :linear,
                  rock_type = nothing, hysteresis = EndpointScaling(), perm_y = nothing,
                  geometry = ve_column_geometry(mesh, columns), tolerance = 1e-6,
                  tip_saturation = 1e-6, touched_saturation = 1e-3) -> LayeredFringe

A capillary fringe in VE columns whose rock varies with depth: the columns
`columns` of the 3D `mesh`, with the fine cells' porosity `poro` and
x-permeability `perm_x` (per cell of `mesh`, as given to
[`ve_domain`](@ref); `perm_y` only checked, as for
[`SharpInterfaceIntegrated`](@ref)). The fine-scale curves, one set for all
cells, or one per rock type (vectors, with `rock_type` the type of each cell
of `mesh`, 1, 2, …):
- `invPc3D`: the fine-scale brine saturation as a function of capillary
  pressure (Pa); a physical curve (for co2lab's dimensionless `(C, alpha)`
  see [`co2lab_capillary_curve`](@ref));
- `kr3D`: the fine-scale CO2 relperm of the CO2 saturation, or an exponent
  `β` for `s^β`;
- `krw3D`: the fine-scale brine relperm of the brine saturation, or
  `:linear` for `(s_w − r_w)/(1 − r_w)`;
- `nothing` for `invPc3D` is a rock type without capillarity (a sharp
  interface in it);
- `hysteresis`: the [`HysteresisLaw`](@ref), [`EndpointScaling`](@ref)
  (co2lab's default) or [`FixedResidual`](@ref);
- `tolerance`, `tip_saturation`, `touched_saturation`: as for
  [`CapillaryFringeTable`](@ref).

Saturations are integrated by pore volume, relperms by permeability, over
the fringe's profile cell by cell, with the hysteresis law's residual CO2.
The caprock's rugosity (the domain's `:rugosity`, `dh` per column) is the
height of its subscale relief: the CO2 there, at the fringe's saturation, is
immobile, and both relperms are integrated from `dh` down (at zero
capillarity co2lab's accretion layer, as [`SharpInterfaceIntegrated`](@ref)).

Returns a new `LayeredFringe` to give to [`co2lab_ve_fluid`](@ref) as
`vertical`; it belongs to this grid. The arguments are only read.
"""
struct LayeredFringe{W, L <: HysteresisLaw} <: VerticalModel
    cells::SharpInterfaceIntegrated          # per column: cell boundaries, ∫φA, ∫k_x dz
    fringes::Vector{<:CapillaryFringeTable}  # per rock type: curves, table settings
    krw3D::Vector{W}                         # per rock type: brine relperm, or :linear
    cell_type::Vector{Int}                   # rock type of cell a of `cells`
    law::L                                   # hysteresis law
end

function LayeredFringe(mesh, columns::VEColumns; poro, perm_x, invPc3D, kr3D = 3,
                       krw3D = :linear, rock_type = nothing,
                       hysteresis::HysteresisLaw = EndpointScaling(), perm_y = nothing,
                       geometry = ve_column_geometry(mesh, columns),
                       tolerance = 1e-6, tip_saturation = 1e-6,
                       touched_saturation = 1e-3)
    types = something(rock_type, ones(Int, number_of_cells(mesh)))
    n = maximum(types)
    curves = [_per_type(x, n) for x in (invPc3D, kr3D, krw3D)]
    for (pc, kw) in zip(curves[1], curves[3])
        pc isa Tuple &&
            throw(ArgumentError("LayeredFringe takes physical capillary curves; " *
                                "co2lab_capillary_curve converts co2lab's (C, alpha)"))
        kw isa Symbol && kw != :linear &&
            throw(ArgumentError("krw3D: a function of the brine saturation, or :linear"))
    end
    cells = SharpInterfaceIntegrated(mesh, columns; poro, perm_x, perm_y, geometry)
    # (nothing: no capillarity, a step at p = 0 to the residual brine)
    fringes = [CapillaryFringeTable(; invPc3D = pc, kr3D = kr, tolerance,
                                    tip_saturation, touched_saturation)
               for (pc, kr) in zip(curves[1], curves[2])]
    krw = [kw for kw in curves[3]]
    return LayeredFringe(cells, fringes, krw, _cell_types(cells, columns, types),
                         hysteresis)
end

# One value per rock type: a vector as given, a single value for all n
_per_type(x, n) = x isa AbstractVector ? (length(x) == n ? x :
                  throw(ArgumentError("one curve per rock type ($n) expected"))) :
                  fill(x, n)

# The rock type of each cell a of `cells` (indexed as its boundaries; the
# last entry of each column unused)
function _cell_types(cells::SharpInterfaceIntegrated, columns::VEColumns, types)
    out = zeros(Int, length(cells.z))
    for i in 1:number_of_columns(columns), (j, c) in enumerate(column_cells(columns, i))
        out[cells.pos[i] + j - 1] = types[c]
    end
    return out
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
columns' rock, a fringe table per rock type, the residual saturations and the
hysteresis law.
"""
struct LayeredFringeModel{T, W, L} <: VerticalModel
    cells::SharpInterfaceIntegrated  # per column: cell boundaries, ∫φ·area, ∫k_x dz
    tables::Vector{T}                # per rock type: integrals over capillary pressure
    krw3D::Vector{W}                 # per rock type: brine relperm, or :linear
    cell_type::Vector{Int}           # rock type of cell a of `cells`
    res_water::Float64               # residual brine saturation
    res_gas::Float64                 # residual CO2 saturation
    law::L                           # hysteresis law
    p_res::Vector{Float64}           # per rock type: capillary pressure where s = r_g
end

# The model belongs to its grid; the tables as CapillaryFringeTable's
function _vertical_model(spec::LayeredFringe, domain, fluid)
    _vertical_model(spec.cells, domain, fluid)
    tables = [_fringe_table(f, domain, fluid) for f in spec.fringes]
    p_res = [_residual_pressure(t, fluid.res_gas) for t in tables]
    return LayeredFringeModel(spec.cells, tables, spec.krw3D, spec.cell_type,
                              fluid.res_water, fluid.res_gas, spec.law, p_res)
end

# The capillary pressure at which table t's drainage CO2 saturation reaches
# rg (0 when the fringe holds rg from its tip on, or rg = 0)
function _residual_pressure(t, rg)
    sg(p) = 1 - t.invpc(max(p, 0.0))
    (rg <= 0 || sg(0.0) >= rg - 1e-12) && return 0.0
    lo, hi = 0.0, t.p[end]
    sg(hi) >= rg ||
        throw(ArgumentError("residual CO2 saturation $rg beyond the capillary curve's"))
    for _ in 1:100
        mid = (lo + hi) / 2
        sg(mid) >= rg ? (hi = mid) : (lo = mid)
    end
    return hi
end

# Cell a's table and brine curve; the settings common to all rock types
_table(m::LayeredFringeModel, a) = m.tables[m.cell_type[a]]
_krw(m::LayeredFringeModel, a) = m.krw3D[m.cell_type[a]]
_settings(m::LayeredFringeModel) = first(m.tables)

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
the free (capillary height `h`), maximum and actual profiles times the pore
volume (m³; actual by the hysteresis law, `_cell_integrals`), `kfree`,
`kmax` the free and maximum profiles' CO2 saturation times `k_x dz`, `krG`,
`krW` the CO2 and brine relperms times `k_x dz` (m³), `pv` and `kz` the
interval's pore volume and `∫ k_x dz`. The column's quantities are those over
`[0, H]`; roadmap Stage 12's transition flux asks for others. Arguments may
carry AD derivatives. Returns a new NamedTuple.
"""
function _interval_integrals(m::LayeredFringeModel, i, T, B, h, h_max, drho)
    c = m.cells
    dg = drho * Jutul.gravity_constant
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
            krW += k * (bot - top) * _brine_relperm_full(_krw(m, a))
            continue
        end

        sf, sm, s, kg = _cell_integrals(m.law, m, a, top, bot, h, h_max, dg)
        free, smax, co2 = free + phiA * sf, smax + phiA * sm, co2 + phiA * s
        kf, km = kf + k * sf, km + k * sm
        krG += k * kg
        krW += k * _brine_integral(m, a, top, bot, s, h, h_max, dg)
    end
    return (free = free, max = smax, co2 = co2, kfree = kf, kmax = km, krG = krG,
            krW = krW, pv = pv, kz = kz)
end

# Cell a's integrals over [top, bot] (within it): ∫ s dz of the free (h) and
# maximum (h_max) drainage profiles and of the actual saturation, and ∫ kr dz
# of the CO2 relperm, by the hysteresis law
function _cell_integrals(::EndpointScaling, m, a, top, bot, h, h_max, dg)
    t, C = _table(m, a), m.res_gas / (1 - m.res_water)
    sf = _depth_integral(t, t.SP, t.s, h, top, bot, dg)
    sm = _depth_integral(t, t.SP, t.s, h_max, top, bot, dg)
    kg = _depth_integral(t, t.krP, t.kr, h, top, bot, dg)
    return (sf, sm, (1 - C) * sf + C * sm, kg)
end

# A fixed residual: the profile of h down to where it holds r_g (depth
# h − d_r, d_r = p_res/(Δρ g)), r_g in the swept zone below down to where the
# old profile held r_g (h_max − d_r), the old profile below (all trapped);
# the CO2 relperm the profile's where it is mobile
function _cell_integrals(::FixedResidual, m, a, top, bot, h, h_max, dg)
    t = _table(m, a)
    clip(x) = _min_c(_max_c(x, top), bot)
    dr = m.p_res[m.cell_type[a]] / dg
    b, e, d = clip(h - dr), clip(h_max - dr), clip(h_max)
    sf = _depth_integral(t, t.SP, t.s, h, top, bot, dg)
    sm = _depth_integral(t, t.SP, t.s, h_max, top, bot, dg)
    s = _depth_integral(t, t.SP, t.s, h, top, b, dg) + m.res_gas * (e - b) +
        _depth_integral(t, t.SP, t.s, h_max, e, d, dg)
    kg = _depth_integral(t, t.krP, t.kr, h, top, b, dg)
    return (sf, sm, s, kg)
end

# The brine relperm of brine alone (s_w = 1), for brine curve `kw`
_brine_relperm_full(kw::Symbol) = 1.0
_brine_relperm_full(kw) = kw(1.0)

# ∫ krw3D(1 − s(z)) dz over [top, bot] of cell a (`s` the actual profile's
# ∫ s dz there): linear, (s_w − r_w)/(1 − r_w), in closed form
function _brine_integral(m::LayeredFringeModel, a, top, bot, s, h, h_max, dg)
    kw = _krw(m, a)
    kw isa Symbol && return (bot - top) - s / (1 - m.res_water)
    return _brine_quadrature(_table(m, a), kw, m, top, bot, h, h_max, dg,
                             m.p_res[m.cell_type[a]])
end

# A function krw3D: Gauss–Legendre between the profiles' kinks (where the
# capillary pressure of either reaches 0 or the entry pressure). The kinks
# carry their AD derivatives: the saturation may jump there (an entry
# pressure, no capillarity), and the moving jump is part of the derivative.
function _brine_quadrature(t, krw3D, m, top, bot, h, h_max, dg, p_res)
    sg(p) = 1 - t.invpc(_max_c(p, 0.0))
    sw(z) = 1 - _local_co2(m.law, m, sg(dg * (h - z)), sg(dg * (h_max - z)))
    z0 = zero(h) + zero(h_max) + zero(dg)
    # where a profile's capillary pressure reaches 0, the entry pressure or (a
    # fixed residual) the pressure holding r_g
    ps = m.law isa FixedResidual ? (0.0, t.entry, p_res) : (0.0, t.entry)
    kinks = [z0 + x - e / dg for x in (h, h_max) for e in ps]
    inside = filter(z -> top < Jutul.value(z) < bot, kinks)
    cuts = sort!([z0 + top; inside; z0 + bot]; by = Jutul.value)
    return sum(_gl_integral(z -> krw3D(sw(z)), a, b)
               for (a, b) in zip(cuts[1:end-1], cuts[2:end]))
end

# The actual CO2 saturation at one depth from the free and maximum profiles'
# there, by the hysteresis law
_local_co2(::EndpointScaling, m, sf, sm) =
    (1 - m.res_gas / (1 - m.res_water)) * sf + m.res_gas / (1 - m.res_water) * sm
_local_co2(::FixedResidual, m, sf, sm) = _max_c(sf, _min_c(m.res_gas, sm))

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
    tip = _settings(m).tip
    Jutul.value(S) < tip && return _layered_newton(m, i, tip, drho) * S / tip
    return _layered_newton(m, i, S, drho)
end

# The root on plain values, then one Newton step with the AD values
function _layered_newton(m, i, S, drho)
    # no buoyancy (CO2 as dense as brine: an unphysical state Newton may pass
    # through, or Jutul's sparsity detection's placeholder state): no height,
    # as CapillaryFringeTable's
    Jutul.value(drho) > 0 || return zero(S) + zero(drho)
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
    c = m.cells
    dg = drho * Jutul.gravity_constant
    v = zero(h) + zero(dg)
    for a in _cell_range(c, i)
        c.z[a] > Jutul.value(h) && break
        dz = c.z[a+1] - c.z[a]
        dz > 0 || continue
        phiA = (c.V[a+1] - c.V[a]) / dz
        t = _table(m, a)
        v += phiA * _depth_integral(t, t.SP, t.s, h, c.z[a], c.z[a+1], dg)
    end
    return v
end

# d/dh of the drainage profile's pore volume of CO2: Σ φA (s(top) − s(bottom))
function _profile_slope(m, i, h, drho)
    c = m.cells
    dg = drho * Jutul.gravity_constant
    slope(t, x) = _table_slope(t, t.SP, t.s, x)
    d = 0.0
    for a in _cell_range(c, i)
        c.z[a] > h && break
        dz = c.z[a+1] - c.z[a]
        dz > 0 || continue
        phiA = (c.V[a+1] - c.V[a]) / dz
        t = _table(m, a)
        d += phiA * (slope(t, dg * (h - c.z[a])) - slope(t, dg * (h - c.z[a+1])))
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
    pend = maximum(t.p[end] for t in m.tables)
    lo, hi = 0.0, pend / (drho * Jutul.gravity_constant) + H
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

# Capillary heights of the mobile and of the maximum profile: endpoint scaling
# the free saturation's; a fixed residual from the column's CO2
_column_heights(m::LayeredFringeModel, i, sG, sGmax, H, drho, rw, rg) =
    _column_heights(m.law, m, i, sG, sGmax, drho, rw, rg)

function _column_heights(::EndpointScaling, m, i, sG, sGmax, drho, rw, rg)
    h = _layered_height(m, i, _free_saturation(sG, sGmax, rw, rg), drho)
    return (h, _layered_height(m, i, sGmax, drho))
end

function _column_heights(::FixedResidual, m, i, sG, sGmax, drho, rw, rg)
    h_max = _layered_height(m, i, sGmax, drho)
    return (_fixed_height(m, i, sG, h_max, drho), h_max)
end

"""
    _fixed_height(m, i, S, h_max, drho) -> Real

A fixed residual's capillary height of column `i` holding the CO2 saturation
`S` after the maximum profile at `h_max`: the root of the column's CO2 less
`S·V`, monotone in the height (its slope `φ·(s_top − r_g)` ≥ 0 wherever the
mobile profile holds more than `r_g`), on plain values to round-off, the
derivatives by one Newton step with the AD values. At or above the maximum's
CO2, `h_max`. Below the CO2 of height 0 (all trapped; dissolution can take
more), linear in `S` from the largest height still holding that much (the
mobile profile then holds none above `r_g`). Zero for `drho ≤ 0`.
"""
function _fixed_height(m::LayeredFringeModel, i, S, h_max, drho)
    z0 = zero(S) + zero(h_max) + zero(drho)
    Jutul.value(drho) > 0 || return z0
    Hc, V = _column_extent(m.cells, i)
    d, hm, v = Jutul.value(drho), Jutul.value(h_max), Jutul.value(S) * V
    F(h) = _interval_integrals(m, i, 0.0, Hc, h, hm, d).co2
    v >= F(hm) && return z0 + h_max
    co2(h) = _interval_integrals(m, i, 0.0, Hc, h, h_max, drho).co2
    F0 = F(0.0)
    if v <= F0
        F0 > 0 || return z0
        return _plateau_end(F, F0, hm) * S * V / co2(zero(hm))
    end
    h = _monotone_root(F, v, 0.0, hm)
    # the slope by central differences within [0, h_max]
    e = 1e-7 * max(hm, 1e-3)
    a, b = max(h - e, 0.0), min(h + e, hm)
    dF = (F(b) - F(a)) / (b - a)
    dF > 0 || return z0 + h
    return h - (co2(h) - S * V) / dF
end

# The largest height at which the column's CO2 F (nondecreasing) is still F0,
# its value at 0, within [0, hm] (plain values; bisection)
function _plateau_end(F, F0, hm)
    tol = 1e-12 * max(F0, eps())
    F(hm) <= F0 + tol && return hm
    lo, hi = 0.0, hm
    for _ in 1:100
        mid = (lo + hi) / 2
        F(mid) <= F0 + tol ? (lo = mid) : (hi = mid)
        hi - lo <= 4eps(hi) && break
    end
    return lo
end

# The root of F(h) = v for F nondecreasing with F(lo) < v < F(hi) (plain
# values): regula falsi (Illinois) safeguarded by bisection, to round-off
function _monotone_root(F, v, lo, hi)
    flo, fhi = F(lo) - v, F(hi) - v
    side = 0
    h = lo
    for _ in 1:200
        h = flo != fhi ? (lo * fhi - hi * flo) / (fhi - flo) : (lo + hi) / 2
        (lo < h < hi) || (h = (lo + hi) / 2)
        r = F(h) - v
        abs(r) <= 1e-13 * max(v, eps()) && return h
        if r > 0
            hi, fhi = h, r
            side == 1 && (flo /= 2)
            side = 1
        else
            lo, flo = h, r
            side == -1 && (fhi /= 2)
            side = -1
        end
        hi - lo <= 4eps(hi) && return h
    end
    return h
end

# Both relperms integrated over the column below the subscale relief (the
# top dh, which takes no part in lateral flow), weighted by permeability and
# scaled by the whole column's ∫ k dz; the endpoints krw, krg are the
# fine-scale curves' (unused)
function _column_relperm(m::LayeredFringeModel, i, sG, sGmax, h, h_max, H, drho, dh,
                         krw, krg, rw, rg)
    Hc, _ = _column_extent(m.cells, i)
    lost = min(dh, Hc)
    I = _interval_integrals(m, i, lost, Hc, h, h_max, drho)
    kz = lost > 0 ? _interval_integrals(m, i, 0.0, Hc, h, h_max, drho).kz : I.kz
    krW = _column_brine(m.law, m, i, sG, sGmax, h_max, drho, lost, Hc, I, rw, rg)
    return (krW / kz, I.krG / kz)
end

# The brine relperm (times ∫ k dz) over [top, Hc] from its integrals I, with
# the tip's treatment (below): endpoint scaling at the free saturation
_column_brine(::EndpointScaling, m, i, sG, sGmax, h_max, drho, top, Hc, I, rw, rg) =
    _tip_brine_relperm(m, i, _free_saturation(sG, sGmax, rw, rg), sGmax, h_max, drho,
                       top, Hc, I)

# A fixed residual: the column below the tip saturation along the ray to
# S_max = tip (as endpoint scaling's whole-column case; the actual
# saturation, not a free one, sets the height)
function _column_brine(::FixedResidual, m, i, sG, sGmax, h_max, drho, top, Hc, I, rw,
                       rg)
    tip, sm = _settings(m).tip, Jutul.value(sGmax)
    sm >= tip && return I.krW
    ratio = sm > 1e-12 * tip ? sG / sGmax : (sm > 0 ? Jutul.value(sG) / sm : 1.0)
    ht = _layered_newton(m, i, tip, drho)
    hr = _fixed_height(m, i, tip * ratio, ht, drho)
    J = _interval_integrals(m, i, top, Hc, hr, ht, drho)
    return I.kz - (I.kz - J.krW) * sGmax / tip
end

# The brine relperm (times ∫ k dz) below the tip saturation. There h is linear
# in S, so the profile at h holds CO2 ∝ S², not S, and the integrated brine
# relperm would have no slope at S = 0. Its deficit (∫ k dz less it) is taken
# linear in each profile's saturation instead, exact at the tip: for the
# linear krw3D the deficit is linear in the profiles' CO2, so this is exact
# (co2lab's krW_simple in a homogeneous column); for another krw3D to first
# order. `I` holds the integrals over the depths [top, Hc] (below the
# subscale relief), and so does the result.
function _tip_brine_relperm(m::LayeredFringeModel, i, Sf, Sm, h_max, drho, top, Hc,
                            I)
    tip = _settings(m).tip
    sf, sm = Jutul.value(Sf), Jutul.value(Sm)
    (sf >= tip && sm >= tip) && return I.krW
    ht = _layered_newton(m, i, tip, drho)
    at(hf, hm) = _interval_integrals(m, i, top, Hc, hf, hm, drho)
    if all(kw -> kw isa Symbol, m.krw3D)
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
# bottoms `zb`; the rock is constant within a cell, so plain depth averages),
# the actual one by the hysteresis law
function _fine_saturations(m::LayeredFringeModel, i, zb, h, h_max, drho, rw, rg)
    dg = drho * Jutul.gravity_constant
    tops = [0.0; zb[1:end-1]]
    out = [b > top ? _cell_integrals(m.law, m, a, top, b, h, h_max, dg) ./ (b - top) :
           (0.0, 0.0, 0.0, 0.0)
           for (a, top, b) in zip(_cell_range(m.cells, i), tops, zb)]
    return ([x[3] for x in out], [x[2] for x in out])
end

# The pore volume below the touched depth (the deepest the maximum profile's
# CO2 reaches; below the touched saturation linear in smax, as
# CapillaryFringeTable's)
function _untouched_brine(m::LayeredFringeModel, i, smax, h_max, H, drho)
    c, touched = m.cells, _settings(m).touched
    if Jutul.value(smax) < touched
        ht = _layered_height(m, i, touched, drho)
        h_t = _touched_depth(m, i, ht, drho) * smax / touched
    else
        h_t = _touched_depth(m, i, h_max, drho)
    end
    Hc, V = _column_extent(c, i)
    return (V - _integral(c, i, c.V, _min_c(h_t, Hc))) / V
end

# The deepest point the drainage profile of capillary height hh reaches:
# where the capillary pressure Δρ g (hh − z) falls to the entry pressure of
# the rock there (cell by cell, top first)
function _touched_depth(m::LayeredFringeModel, i, hh, drho)
    c = m.cells
    dg = drho * Jutul.gravity_constant
    d = zero(hh) + zero(dg)
    for a in _cell_range(c, i)
        reach = hh - _table(m, a).entry / dg
        Jutul.value(reach) > c.z[a] || continue
        d = _min_c(reach, c.z[a+1])
    end
    return _max_c(d, 0.0)
end

# The inventory's split as CapillaryFringeTable's, each profile's fraction
# between two depths by pore volume; co2lab_compatible: co2lab's by height
function _category_saturations(m::LayeredFringeModel, S, Smax, h, h_max, col, drho,
                               rw, rg; co2lab_compatible::Bool = false)
    co2lab_compatible && return _height_split(S, Smax, h, h_max, col, rw, rg)
    m.law isa FixedResidual && return _fixed_categories(m, S, h, h_max, col, drho)
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

# The inventory's split with a fixed residual: the immobile part is the
# residual min(r_g, s_max) the column keeps when its mobile CO2 is gone (the
# profile at h = 0), the mobile part the rest, (s_h − r_g)⁺ at each depth;
# each part's share of a zone from its own depth distribution. "Below the
# plume" is below the capillary height h, as for endpoint scaling.
function _fixed_categories(m::LayeredFringeModel, S, h, h_max, col, drho)
    n = eachindex(S)
    I(i, hh, a, b) = (Hc = _column_extent(m.cells, i)[1];
                      _interval_integrals(m, i, clamp(a, 0.0, Hc), clamp(b, 0.0, Hc), hh,
                                          h_max[i], drho[i]).co2)
    V(i) = _column_extent(m.cells, i)[2]
    # each part's pore volume of CO2 between depths a and b
    imm(i, a, b) = I(i, 0.0, a, b)
    mob(i, a, b) = I(i, h[i], a, b) - I(i, 0.0, a, b)
    S_immob = [min(imm(i, 0.0, col.H[i]) / V(i), S[i]) for i in n]
    S_mob = S .- S_immob
    # a part's share of the zone [a, b] (0 for an empty part)
    share(f, i, a, b) = (t = f(i, 0.0, col.H[i]); t > 0 ? max(f(i, a, b), 0.0) / t : 0.0)
    dh, zt = col.dh, col.zt
    struct_end(i) = max(zt[i], dh[i])
    immob_sub = [S_immob[i] * share(imm, i, 0.0, dh[i]) for i in n]
    mobile_sub = [S_mob[i] * share(mob, i, 0.0, dh[i]) for i in n]
    immob_struct = [S_immob[i] * share(imm, i, dh[i], struct_end(i)) for i in n]
    mobile_struct = [S_mob[i] * share(mob, i, dh[i], struct_end(i)) for i in n]
    res = [S_immob[i] * share(imm, i, max(h[i], zt[i]), col.H[i]) for i in n]
    return (immob_sub = immob_sub, mobile_sub = mobile_sub,
            immob_struct = immob_struct, mobile_struct = mobile_struct,
            immob_free = S_immob .- immob_struct .- immob_sub,
            mobile_free = S_mob .- mobile_struct .- mobile_sub, res = res)
end

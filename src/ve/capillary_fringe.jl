# co2lab's capillary-fringe vertical models, 'P-scaled table' and (corrected)
# 'P-K-scaled table' (co2lab: addVERelpermCapillaryFringe, makeVEFluid).
#
# Below the caprock the fine-scale capillary pressure is linear in depth,
# pc(z) = Pc − Δρ g z, so a column's integrals become integrals over capillary
# pressure. One table, built once per model, holds them as functions of p:
#   SP(p)  = ∫₀^p sg(u) du,    krP(p) = ∫₀^p kr3D(sg(u)) du,
# with sg(u) = 1 − invPc3D(u) the fine-scale CO2 saturation. A column with CO2
# saturation S, height H and density difference Δρ then has
#   S · H · Δρ g = SP(Pc) − SP(Pc − Δρ g H)     (2nd term: plume at the bottom),
#   krG = (krP(Pc) − krP(Pc − Δρ g H)) / (H Δρ g),   capillary height h = Pc/(Δρ g).
#
# The table: cubic Hermite on adaptively placed nodes with exact values and
# derivatives (a node at the entry pressure). About 250 nodes at tolerance 1e-6
# beat co2lab's 200 000 uniform, linearly interpolated samples, which are poor
# for thin plumes. co2lab's range (Pmax) and its scaling of a (C, alpha) curve
# are kept.
#
# The plume tip: dPc/dS is unbounded as S leaves 0 (Pc grows like √S, or jumps
# to the entry pressure), so below a small free saturation Pc is linear in S.
#
# Leverett scaling ('P-K-scaled table'): the curve is a J-curve of
# J = pc/kscale, kscale = σ·√(φ/k) per column. The table is over J and a column
# works in J with the scaled density difference Δρ/kscale, so inversion,
# truncation and krG stay consistent and h comes out physical. co2lab's own
# version fails on first use and mixes J with pressures; this one is corrected.

# ---------------------------------- Settings ----------------------------------

"""
    CapillaryFringeTable(; invPc3D, kr3D = 3, tolerance = 1e-6,
                           tip_saturation = 1e-6, touched_saturation = 1e-3,
                           leverett = false, surface_tension = 30e-3)

co2lab's `'P-scaled table'` vertical model: a capillary fringe between plume
and brine, with one fine-scale capillary curve for the whole formation. Give
it to [`co2lab_ve_fluid`](@ref) as `vertical`; it is grid-independent, and
[`setup_co2lab_model`](@ref) builds the table for the model's domain.

# Fields

Each is also a constructor keyword of the same name, and makeVEFluid's option
of that name.

- `invPc3D`: fine-scale brine saturation as a function of capillary pressure
  (Pa); or `(C, alpha)` for co2lab's
  `sw = max((C'/(p + C'))^(1/alpha), res_water)` with `C' = C·g·Hmax·Δρ_ref`,
  scaled as in co2lab by the grid's tallest column `Hmax` and the reference
  density difference.
- `kr3D`: fine-scale CO2 relperm as a function of CO2 saturation, or an
  exponent `β` for `s^β` (co2lab's default 3).
- `tolerance`: the table's relative accuracy. co2lab's uniform table of
  200 000 samples is not offered.
- `tip_saturation`: below this free CO2 saturation Pc is linear in it (at the
  plume tip the exact Pc has an unbounded derivative).
- `touched_saturation`: below this maximum CO2 saturation the depth of brine
  touched by CO2 (dissolution) is linear in it. The exact depth grows like √S
  at the tip, far faster than the free CO2, and Newton cycles there.
- `leverett`: co2lab's `'P-K-scaled table'`, corrected: `invPc3D` is a
  function of `J = pc/kscale`, `kscale = surface_tension·√(φ/k)` per column
  from the domain's upscaled porosity and permeability.
- `surface_tension`: `σ` of the Leverett scaling (N/m).

# Notes

- Does not combine with caprock rugosity (nor does co2lab's):
  [`setup_co2lab_model`](@ref) refuses a domain with non-zero `:rugosity`.
- Hysteresis is co2lab's `'endpoint scaling'` (its default): the residual
  CO2 at each depth scales with the maximum saturation reached there.
  co2lab's experimental `'fixed residual'` is not offered.
"""
struct CapillaryFringeTable{F, K} <: VerticalModel
    invPc3D::F                   # brine saturation of capillary pressure, or (C, alpha)
    kr3D::K                      # CO2 relperm of CO2 saturation, or a Corey exponent
    tolerance::Float64           # relative accuracy of the table
    tip_saturation::Float64      # free saturation below which Pc is linear in it
    touched_saturation::Float64  # max. saturation below which touched depth is linear
    leverett::Bool               # Leverett scaling ('P-K-scaled table')
    surface_tension::Float64     # σ of the Leverett scaling (N/m)
end

CapillaryFringeTable(; invPc3D, kr3D = 3, tolerance = 1e-6, tip_saturation = 1e-6,
                     touched_saturation = 1e-3, leverett::Bool = false,
                     surface_tension = 30e-3) =
    CapillaryFringeTable(invPc3D, kr3D, Float64(tolerance), Float64(tip_saturation),
                         Float64(touched_saturation), leverett, Float64(surface_tension))

"""
    PScaledTable

The built form of [`CapillaryFringeTable`](@ref) for one domain: the column
integrals over capillary pressure as a cubic Hermite table, shared by all
columns (over `J`, with each column's `kscale`, under Leverett scaling).

# Fields

- `p`: the table's nodes, capillary pressure (Pa), or `J` under Leverett
  scaling, from 0 to the upper limit.
- `SP`: integral of the fine-scale CO2 saturation over `p`, at each node.
- `s`: its derivative, the fine-scale CO2 saturation at each node.
- `krP`: integral of the fine-scale CO2 relperm over `p`, at each node.
- `kr`: its derivative, the fine-scale CO2 relperm at each node.
- `invpc`: fine-scale brine saturation as a function of capillary pressure
  (or of `J`).
- `kscale`: per column, the Leverett factor `σ√(φ/k)` (Pa); 1 without
  Leverett scaling.
- `tip`: tip saturation: below this free CO2 saturation, Pc is linear in it.
- `entry`: entry pressure (or `J`): the fine-scale CO2 saturation is positive
  above it; 0 without entry pressure.
- `touched`: maximum CO2 saturation below which the depth of brine touched
  by CO2 is linear in it.
"""
struct PScaledTable{F} <: VerticalModel
    p::Vector{Float64}       # nodes (capillary pressure, Pa; or J), 0 to Pmax
    SP::Vector{Float64}      # integral of the fine-scale CO2 saturation over p
    s::Vector{Float64}       # its derivative: the fine-scale CO2 saturation
    krP::Vector{Float64}     # integral of the fine-scale CO2 relperm over p
    kr::Vector{Float64}      # its derivative: the fine-scale CO2 relperm
    invpc::F                 # fine-scale brine saturation of capillary pressure (or J)
    kscale::Vector{Float64}  # Leverett factor σ√(φ/k) per column (Pa); 1 without
    tip::Float64             # tip saturation S_r
    entry::Float64           # entry pressure (or J): sg > 0 above it; 0 without
    touched::Float64         # max. saturation below which the touched depth is linear
end

# ------------------------------- Building tables ------------------------------

# co2lab's fringe for `domain` and `fluid`: its table; refuses caprock rugosity
function _vertical_model(spec::CapillaryFringeTable, domain, fluid)
    _refuse_rugosity(domain, "CapillaryFringeTable")
    return _fringe_table(spec, domain, fluid)
end

# The table for `domain` and `fluid`, with co2lab's fine-scale functions
# (setup_fine_scale_functions) and upper limit. Returns a new PScaledTable.
function _fringe_table(spec::CapillaryFringeTable, domain, fluid)
    H = domain[:column_height]
    g = Jutul.gravity_constant
    drho_ref = fluid.rho_ref[1] - fluid.rho_ref[2]
    invpc = _fine_capillary_curve(spec.invPc3D, maximum(H) * drho_ref * g,
                                  fluid.res_water)
    kr3D = spec.kr3D isa Number ? (s -> s^spec.kr3D) : spec.kr3D

    # Leverett: each column's factor from its upscaled rock (co2lab's rock2D)
    kscale = spec.leverett ?
             sqrt.(domain[:porosity] ./ domain[:permeability]) .* spec.surface_tension :
             ones(length(H))
    # the upper limit in the table's variable (pressure, or J)
    pmax = _upper_capillary_pressure(invpc, fluid.res_water,
                                     maximum(H .* drho_ref .* g ./ kscale))
    sg(p) = 1 - invpc(max(p, 0.0))
    t = _adaptive_table(sg, kr3D, pmax, spec.tolerance)
    return PScaledTable(t..., invpc, kscale, spec.tip_saturation,
                        _entry_pressure(sg, pmax), spec.touched_saturation)
end

# Throws if the domain has caprock rugosity: the fringe models (like co2lab's)
# ignore it in the flow, while the trapping inventory would split off
# subscale-trapped CO2.
function _refuse_rugosity(domain, model)
    haskey(domain, :rugosity) && any(>(0), domain[:rugosity]) &&
        throw(ArgumentError("$model does not combine with caprock rugosity; " *
                            "set the domain's rugosity to zero"))
    return nothing
end

# The density difference in the table's variable: Δρ/kscale (Δρ without
# Leverett scaling), so that column integrals over p are those over J
_table_drho(m::PScaledTable, i, drho) = drho / m.kscale[i]

# Fraction of column i's brine untouched by CO2. Brine is touched down to
# where the maximum profile's CO2 saturation is positive: h_max less the entry
# height (h_max alone would let a trace of CO2 touch metres of brine). Below
# the touched saturation the depth is linear in smax (the exact one grows like
# √smax, outgrowing the free CO2, and Newton cycles).
function _untouched_brine(m::PScaledTable, i, smax, h_max, H, drho)
    drho_g = _table_drho(m, i, drho) * Jutul.gravity_constant
    if _plain(smax) < m.touched
        pc_t = _pc_newton(m, m.touched, H, _table_drho(m, i, drho))
        h_t = _max_c(pc_t - m.entry, 0.0) / drho_g * smax / m.touched
    else
        h_t = _max_c(h_max - m.entry / drho_g, 0.0)
    end
    return (H - _min_c(h_t, H)) / H
end

# A function, or co2lab's (C, alpha), scaled by the pressure drop over the
# tallest column (`dp_max` = Hmax·Δρ_ref·g)
_fine_capillary_curve(f, dp_max, rw) = f
# no capillarity: a step at p = 0, from brine to the residual brine
_fine_capillary_curve(::Nothing, dp_max, rw) = p -> p > 0 ? rw : 1.0
function _fine_capillary_curve((C, a)::Tuple{Real, Real}, dp_max, rw)
    Cs = C * dp_max
    return p -> max((Cs / (p + Cs))^(1 / a), rw)
end

# co2lab's determine_upper_cap_press_limit: the pressure doubled from 1 Pa
# until the brine is within 1e-3 of residual, plus the drop over the tallest
# column (`dp_max`), so that Pmax saturates every column
function _upper_capillary_pressure(invpc, rw, dp_max)
    dp = 1.0
    while invpc(dp) > rw + 1e-3
        dp *= 2
    end
    return dp + dp_max
end

# 5-point Gauss–Legendre on [-1, 1]
const _GL_X = (-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683,
               0.906179845938664)
const _GL_W = (0.236926885056189, 0.478628670499366, 0.568888888888889,
               0.478628670499366, 0.236926885056189)

# ∫ f over [a, b]: 5-point Gauss–Legendre on 4 sub-intervals
function _gl_integral(f, a, b)
    h = (b - a) / 4
    return sum(h / 2 * sum(w * f(a + (k - 0.5) * h + h / 2 * x)
                           for (x, w) in zip(_GL_X, _GL_W)) for k in 1:4)
end

"""
    _adaptive_table(sg, kr3D, pmax, rtol) -> (p, SP, s, krP, kr)

The column integrals `SP = ∫ sg dp`, `krP = ∫ kr3D(sg) dp` on nodes from 0 to
`pmax`: 17 uniform nodes and the entry pressure (where `sg` turns positive),
each interval halved until the cubic Hermite midpoint values of both integrals
match Gauss–Legendre to `rtol` (relative). Returns new vectors: the nodes, the
integrals and their derivatives at the nodes.
"""
function _adaptive_table(sg, kr3D, pmax, rtol)
    krs(p) = kr3D(sg(p))
    nodes = collect(range(0.0, pmax; length = 17))
    pe = _entry_pressure(sg, pmax)
    pe > 0 && push!(nodes, pe)
    sort!(unique!(nodes))
    todo = collect(zip(nodes[1:end-1], nodes[2:end]))
    done = Tuple{Float64, Float64}[]
    while !isempty(todo)
        a, b = pop!(todo)
        if _hermite_fine(sg, krs, a, b, rtol) || b - a < 1e-9 * pmax
            push!(done, (a, b))
        else
            push!(todo, (a, (a + b) / 2), ((a + b) / 2, b))
        end
    end
    sort!(done)
    p = [first.(done); last(done[end])]
    SP = cumsum([0.0; [_gl_integral(sg, a, b) for (a, b) in done]])
    krP = cumsum([0.0; [_gl_integral(krs, a, b) for (a, b) in done]])
    return (p, SP, _right_values(sg, p), krP, _right_values(krs, p))
end

# f at the nodes p, at 0 from the right: a curve without capillarity steps up
# at p = 0, and the integrals' slope there is the plume's (continuous curves:
# the same value)
_right_values(f, p) = [x == 0 ? f(nextfloat(0.0)) : f(x) for x in p]

# Where sg turns positive (0 without an entry pressure), by bisection
function _entry_pressure(sg, pmax)
    sg(1e-9 * pmax) > 0 && return 0.0
    lo, hi = 0.0, pmax
    for _ in 1:200
        mid = (lo + hi) / 2
        sg(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    return hi
end

# Whether the Hermite interpolant of both integrals on [a, b] matches
# Gauss–Legendre at the midpoint to rtol (relative, with a tiny floor)
function _hermite_fine(sg, krs, a, b, rtol)
    m = (a + b) / 2
    for f in (sg, krs)
        I_ab, I_am = _gl_integral(f, a, b), _gl_integral(f, a, m)
        # the node values as the table has them (from the right at 0)
        fa = a == 0 ? f(nextfloat(0.0)) : f(a)
        pred = _hermite(a, b, 0.0, I_ab, fa, f(b), m)
        abs(pred - I_am) <= rtol * abs(I_am) + 1e-14 * (b - a) || return false
    end
    return true
end

# ---------------------------------- Lookups -----------------------------------

"Cubic Hermite on [x0, x1] with values y0, y1 and derivatives d0, d1."
function _hermite(x0, x1, y0, y1, d0, d1, x)
    dx = x1 - x0
    t = (x - x0) / dx
    return (2t^3 - 3t^2 + 1) * y0 + (t^3 - 2t^2 + t) * dx * d0 +
           (-2t^3 + 3t^2) * y1 + (t^3 - t^2) * dx * d1
end

"Derivative of `_hermite` at x."
function _hermite_slope(x0, x1, y0, y1, d0, d1, x)
    dx = x1 - x0
    t = (x - x0) / dx
    return (6t^2 - 6t) / dx * y0 + (3t^2 - 4t + 1) * d0 +
           (-6t^2 + 6t) / dx * y1 + (3t^2 - 2t) * d1
end

"""
    _table_value(m, Y, dY, x) -> Real

The table `Y` (`m.SP` or `m.krP`, derivatives `dY`) at `x`: 0 for `x < 0`,
cubic Hermite on the nodes (at `x = 0` with the slope from the right, the
plume's), linear beyond the last node. `x` may carry AD derivatives.
"""
function _table_value(m::PScaledTable, Y, dY, x)
    xv = _plain(x)
    xv < 0 && return zero(x)
    n = length(m.p)
    xv >= m.p[n] && return Y[n] + dY[n] * (x - m.p[n])
    k = searchsortedlast(m.p, xv)
    return _hermite(m.p[k], m.p[k+1], Y[k], Y[k+1], dY[k], dY[k+1], x)
end

"The derivative of `_table_value` at the plain value `x`."
function _table_slope(m::PScaledTable, Y, dY, x)
    x < 0 && return 0.0
    n = length(m.p)
    x >= m.p[n] && return dY[n]
    k = searchsortedlast(m.p, x)
    return _hermite_slope(m.p[k], m.p[k+1], Y[k], Y[k+1], dY[k], dY[k+1], x)
end

_SP(m::PScaledTable, x) = _table_value(m, m.SP, m.s, x)
_krP(m::PScaledTable, x) = _table_value(m, m.krP, m.kr, x)

# ------------------------------- Column integrals -----------------------------

"""
    _free_saturation(sg, sGmax, rw, rg) -> Real

co2lab's `free_sg`: the part of `sg` that is mobile and will not stay behind
as residual CO2 (endpoint scaling; Nilsen et al. 2016, eq. 17).
"""
function _free_saturation(sg, sGmax, rw, rg)
    fsg = ((1 - rw) * sg - sGmax * rg) / (1 - rw - rg)
    # drainage (no residual trapping): the saturation itself
    _plain(sg) >= _plain(sGmax) && (fsg = sg)
    return _max_c(fsg, 0.0)
end

"""
    _fringe_capillary_pressure(m, sg, H, drho) -> Pc

The upscaled capillary pressure of a column holding the free CO2 saturation
`sg` (density difference `drho` in the table's variable): the root of
`SP(Pc) − SP(Pc − Δρ g H) − sg·H·Δρ g`, the second term only for a plume
reaching the bottom. Solved on plain values to round-off; the derivatives by
one Newton step with the AD values. Below the tip saturation, linear in `sg`
from 0 to its value there.
"""
function _fringe_capillary_pressure(m::PScaledTable, sg, H, drho)
    sv = _plain(sg)
    if sv < m.tip
        # the plume tip: linear in the saturation (Pc at the tip keeps Δρ's
        # derivatives)
        return _pc_newton(m, m.tip, H, drho) * sg / m.tip
    end
    return _pc_newton(m, sg, H, drho)
end

# The root on plain values, then one Newton step with the AD values
function _pc_newton(m, sg, H, drho)
    dP = _plain(drho) * Jutul.gravity_constant * H
    pc = _solve_pc(m, _plain(sg) * dP, dP)
    dPad = drho * Jutul.gravity_constant * H
    F = _SP(m, pc) - _SP(m, pc - dPad) - sg * dPad
    dF = _table_slope(m, m.SP, m.s, pc) - _table_slope(m, m.SP, m.s, pc - dP)
    # a full column: Pc is not determined (any larger one fills it too), so
    # no derivatives
    dF > 0 || return zero(F) + pc
    return pc - F / dF
end

# Pc with SP(Pc) − SP(Pc − dP) = SP_col (plain values), the plume reaching
# the bottom or not; Newton safeguarded by bisection, to round-off
function _solve_pc(m, SP_col, dP)
    SP_col <= 0 && return 0.0
    f(pc) = _SP(m, pc) - _SP(m, pc - dP) - SP_col
    df(pc) = _table_slope(m, m.SP, m.s, pc) - _table_slope(m, m.SP, m.s, pc - dP)
    # f is nondecreasing and constant beyond p_end + dP (the table saturated
    # over the whole column): more CO2 than the table can hold gives that Pc
    lo, hi = 0.0, m.p[end] + dP
    f(hi) <= 0 && return hi
    pc = (lo + hi) / 2
    for _ in 1:200
        r = f(pc)
        r > 0 ? (hi = pc) : (lo = pc)
        abs(r) <= 1e-13 * max(SP_col, eps()) && return pc
        d = df(pc)
        next = d > 0 ? pc - r / d : (lo + hi) / 2
        pc = lo < next < hi ? next : (lo + hi) / 2
        hi - lo <= 4eps(hi) && return pc
    end
    return pc
end

# ------------------------- Vertical-model interface ---------------------------

# Capillary heights Pc/(Δρ g) of the free saturation and of the maximum
# saturation (co2lab's upscaledSat2height, capillary case)
function _column_heights(m::PScaledTable, i, sG, sGmax, H, drho, rw, rg)
    drho = _table_drho(m, i, drho)
    drho_g = drho * Jutul.gravity_constant
    pc = _fringe_capillary_pressure(m, _free_saturation(sG, sGmax, rw, rg), H, drho)
    pc_max = _fringe_capillary_pressure(m, sGmax, H, drho)
    return (pc / drho_g, pc_max / drho_g)
end

# co2lab's krG_ptable and krW_simple (the fringe's brine relperm, a different
# function from the sharp interface's krW_simple: linear in the brine not
# bound as residual by the free CO2). No accretion layer: co2lab does not
# combine the fringe with rugosity.
function _column_relperm(m::PScaledTable, i, sG, sGmax, h, h_max, H, drho, dh,
                         krw, krg, rw, rg)
    drho_g = _table_drho(m, i, drho) * Jutul.gravity_constant
    pc = h * drho_g
    # the part below the column's bottom (krP is 0 below 0)
    krG = (_krP(m, pc) - _krP(m, pc - drho_g * H)) / (H * drho_g)

    sw = 1 - sG
    sw_eff = sw - _free_saturation(sG, sGmax, rw, rg) / (1 - rw) * rw
    return (_max_c(sw_eff, 0.0), krG)
end

# co2lab's cap_fringe_h2s, endpoint scaling: each fine
# cell's average of the current (free) and the maximum fringe profile, mixed
# as s = (1 − C)·s_free + C·s_max with C = rg/(1 − rw)
function _fine_saturations(m::PScaledTable, i, zb, h, h_max, drho, rw, rg)
    drho = _table_drho(m, i, drho)
    tops = [0.0; zb[1:end-1]]
    s_free = [_cell_fringe_saturation(m, a, b, h, drho) for (a, b) in zip(tops, zb)]
    s_max = [_cell_fringe_saturation(m, a, b, h_max, drho) for (a, b) in zip(tops, zb)]
    C = rg / (1 - rw)
    return ((1 - C) .* s_free .+ C .* s_max, s_max)
end

# The fringe saturation averaged over a cell from depth a to b below the
# caprock, for capillary height h: co2lab's trapezoid rule on 10 levels
function _cell_fringe_saturation(m::PScaledTable, a, b, h, drho)
    n, g = 10, Jutul.gravity_constant
    s = 0.0
    for k in 1:n
        w = k == 1 || k == n ? 0.5 : 1.0
        d = h - (a * (n - k) + b * (k - 1)) / (n - 1)
        # no CO2 below the capillary height
        d > 0 && (s += w * (1 - m.invpc(d * drho * g)))
    end
    return s / (n - 1)
end

# The inventory's split of a fringe column: the mobile part by the current
# saturation profile, the immobile part by the maximum one; co2lab_compatible:
# co2lab's split by height
function _category_saturations(m::PScaledTable, S, Smax, h, h_max, col, drho, rw, rg;
                               co2lab_compatible::Bool = false)
    co2lab_compatible && return _height_split(S, Smax, h, h_max, col, rw, rg)
    S_immob = Smax .* rg ./ (1 - rw)
    S_mob = S .- S_immob
    g = Jutul.gravity_constant
    # in the table's variable (J with Leverett scaling)
    drho = [_table_drho(m, i, drho[i]) for i in eachindex(drho)]
    pc, pc_max = h .* drho .* g, h_max .* drho .* g

    # the fraction of each profile between two depths below the caprock
    mob(i, a, b) = _profile_fraction(m, pc[i], drho[i] * g, col.H[i], a, b)
    imm(i, a, b) = _profile_fraction(m, pc_max[i], drho[i] * g, col.H[i], a, b)
    n = eachindex(S)
    dh, zt = col.dh, col.zt
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

"""
    _profile_fraction(m, pc, drho_g, H, a, b) -> Float64

The fraction of a fringe column's CO2 (capillary pressure `pc` at the caprock,
column height `H`) that lies between depths `a` and `b` below the caprock:
`(SP(pc − Δρ g a) − SP(pc − Δρ g b)) / (SP(pc) − SP(pc − Δρ g H))`, depths and
pressures clipped to the column and to zero. Zero for an empty column.
"""
function _profile_fraction(m::PScaledTable, pc, drho_g, H, a, b)
    SPat(z) = _SP(m, pc - drho_g * clamp(z, 0.0, H))
    total = SPat(0.0) - SPat(H)
    total > 0 || return 0.0
    return max(SPat(a) - SPat(b), 0.0) / total
end

# co2lab's capillary-fringe VE fluid, 'P-scaled table' (roadmap Stage 6 item 1;
# co2lab: addVERelpermCapillaryFringe.m, makeVEFluid.m).
#
# The fine-scale capillary pressure is linear in depth below the caprock,
# pc(z) = Pc − Δρ g z, so the column's integrals become integrals over
# capillary pressure. One table, built once, holds them as functions of p:
#   SP(p)  = ∫₀^p sg(u) du,    krP(p) = ∫₀^p kr3D(sg(u)) du,
# with sg(u) = 1 − invPc3D(u) the fine-scale CO2 saturation. A column with CO2
# saturation S, height H and density difference Δρ then has
#   S · H · Δρ g = SP(Pc) − SP(Pc − Δρ g H),
# the second term only when the plume reaches the bottom (truncated). The
# upscaled capillary pressure Pc is found from S by the inverse table, with a
# Newton solve for truncated plumes; krG = (krP(Pc) − krP(Pc − Δρ g H)) /
# (H Δρ g); the capillary height h = Pc/(Δρ g).
#
# Everything is built as co2lab builds it, so that the gates compare like
# with like: the sample count, the grid-dependent upper limit Pmax, the
# scaling of a numeric (C, alpha) curve by the tallest column, the smoothing
# of an entry pressure, the midpoint cumulative sum and linear interpolation
# with extrapolation (MRST's interpTable; at a node the slope of the interval
# to the right).

# ---------------------------------- Settings ----------------------------------

"""
    CapillaryFringeTable(; invPc3D, kr3D = 3, hysteresis = :endpoint_scaling,
                           samples = 200_000)

co2lab's `'P-scaled table'` vertical model: a capillary fringe between the
plume and the brine, with one fine-scale capillary curve for the whole
formation (makeVEFluid's options of the same names):
- `invPc3D`: the fine-scale brine saturation as a function of capillary
  pressure (Pa), or `(C, alpha)` for co2lab's
  `sw = max((C'/(p + C'))^(1/alpha), res_water)`, where `C' = C·g·Hmax·Δρ_ref`
  is scaled by the tallest column `Hmax` of the grid and the reference
  density difference, as co2lab does;
- `kr3D`: the fine-scale CO2 relperm as a function of CO2 saturation, or an
  exponent `β` for `s^β` (co2lab's default 3);
- `hysteresis`: `:endpoint_scaling` (co2lab's default) or `:fixed_residual`;
- `samples`: table samples (makeVEFluid asks for 200 000).

A grid-independent specification; [`setup_co2lab_model`](@ref) builds the
table for the model's domain.
"""
struct CapillaryFringeTable{F, K} <: VerticalModel
    invPc3D::F          # brine saturation of capillary pressure, or (C, alpha)
    kr3D::K             # CO2 relperm of CO2 saturation, or a Corey exponent
    hysteresis::Symbol  # :endpoint_scaling or :fixed_residual
    samples::Int        # table samples
end

function CapillaryFringeTable(; invPc3D, kr3D = 3, hysteresis = :endpoint_scaling,
                              samples = 200_000)
    hysteresis in (:endpoint_scaling, :fixed_residual) ||
        throw(ArgumentError("hysteresis must be :endpoint_scaling or :fixed_residual"))
    return CapillaryFringeTable(invPc3D, kr3D, hysteresis, Int(samples))
end

"""
    PScaledTable

The built form of [`CapillaryFringeTable`](@ref): co2lab's table over
capillary pressure (`make_CO2_table_p_based`), shared by all columns.
"""
struct PScaledTable <: VerticalModel
    p::Vector{Float64}    # capillary pressure axis (Pa), uniform from 0 to Pmax
    SP::Vector{Float64}   # integral of the fine-scale CO2 saturation over p
    krP::Vector{Float64}  # integral of the fine-scale CO2 relperm over p
    hysteresis::Symbol    # :endpoint_scaling or :fixed_residual
end

# ------------------------------- Building tables ------------------------------

# The table for `domain` and `fluid`: co2lab's fine-scale functions
# (setup_fine_scale_functions), upper limit and table
function _vertical_model(spec::CapillaryFringeTable, domain, fluid)
    H = domain[:column_height]
    g = Jutul.gravity_constant
    drho_ref = fluid.rho_ref[1] - fluid.rho_ref[2]
    invpc = _fine_capillary_curve(spec.invPc3D, maximum(H) * drho_ref * g,
                                  fluid.res_water)
    kr3D = spec.kr3D isa Number ? (s -> s^spec.kr3D) : spec.kr3D

    pmax = _upper_capillary_pressure(invpc, fluid.res_water, maximum(H) * drho_ref * g)
    p, SP, krP = _p_based_table(invpc, kr3D, spec.samples, pmax)
    return PScaledTable(p, SP, krP, spec.hysteresis)
end

# A function, or co2lab's (C, alpha), scaled by the pressure drop over the
# tallest column (`dp_max` = Hmax·Δρ_ref·g)
_fine_capillary_curve(f, dp_max, rw) = f
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

# co2lab's make_CO2_table_p_based: samples on a uniform axis, an entry-pressure
# plateau of zero saturation smoothed into a ramp (to keep SP invertible),
# and the midpoint cumulative sums. Returns new vectors (p, SP, krP).
function _p_based_table(invpc, kr3D, samples, pmax)
    p = collect(range(0.0, pmax; length = samples))
    dp = p[2] - p[1]
    sg = 1 .- invpc.(p)
    kr = kr3D.(sg)

    # the zero-saturation plateau below the entry pressure becomes a ramp up
    # to a tenth of the first nonzero value (co2lab's SMOOTH_FAC)
    last_zero = findlast(==(0.0), sg)
    if !isnothing(last_zero) && last_zero > 1
        sg[1:last_zero] = range(0.0, 0.1 * sg[last_zero + 1]; length = last_zero)
    end
    SP = (cumsum(sg) .- sg ./ 2) .* dp
    krP = (cumsum(kr) .- kr ./ 2) .* dp
    return (p, SP, krP)
end

# ---------------------------------- Lookups -----------------------------------

"""
    _interp_table(X, Y, x) -> Real

MRST's `interpTable`: linear interpolation in the increasing table `X → Y`,
extrapolated linearly beyond its ends; at a node, the slope of the interval
to the right (`dinterpq1`), so that AD derivatives match MRST's.
"""
function _interp_table(X, Y, x)
    n = length(X)
    k = clamp(searchsortedlast(X, Jutul.value(x)), 1, n - 1)
    t = (x - X[k]) / (X[k+1] - X[k])
    return Y[k] + t * (Y[k+1] - Y[k])
end

# Slope of the table X → Y at x (the derivative `_interp_table` carries)
function _interp_slope(X, Y, x)
    n = length(X)
    k = clamp(searchsortedlast(X, x), 1, n - 1)
    return (Y[k+1] - Y[k]) / (X[k+1] - X[k])
end

# ------------------------------- Column integrals -----------------------------

"""
    _free_saturation(sg, sGmax, rw, rg) -> Real

co2lab's `free_sg`: the part of `sg` that is mobile and will not stay behind
as residual CO2 (endpoint scaling; Nilsen et al. 2016, eq. 17).
"""
function _free_saturation(sg, sGmax, rw, rg)
    fsg = ((1 - rw) * sg - sGmax * rg) / (1 - rw - rg)
    # drainage (no residual trapping): the saturation itself
    sg >= sGmax && (fsg = sg)
    return _max_c(fsg, 0.0)
end

"""
    _fringe_capillary_pressure(m, sg, H, drho) -> (Pc, pb)

The upscaled capillary pressure `Pc` of a column holding the free CO2
saturation `sg` (co2lab's `compute_untruncated_SP` and `pcWG_ptable`), and
the capillary pressure `pb` at the column's bottom when the plume reaches it
(0 otherwise). A truncated plume holds `SP(Pc) − SP(pb)` with
`Pc = pb + Δρ g H`; `pb` is found by Newton on plain values, its derivatives
by one Newton step with the AD values.
"""
function _fringe_capillary_pressure(m::PScaledTable, sg, H, drho)
    drho_g = drho * Jutul.gravity_constant
    SP = sg * H * drho_g
    pc = _interp_table(m.SP, m.p, SP)
    # the plume fits in the column: no truncation
    pc <= drho_g * H && return (pc, zero(pc))

    pb = _bottom_capillary_pressure(m, Jutul.value(SP), Jutul.value(drho_g * H),
                                    Jutul.value(pc - drho_g * H))
    # one Newton step with AD values: the derivatives of the implicit pb
    f = _interp_table(m.p, m.SP, pb) + SP - _interp_table(m.p, m.SP, pb + drho_g * H)
    dP = Jutul.value(drho_g * H)
    df = _interp_slope(m.p, m.SP, pb) - _interp_slope(m.p, m.SP, pb + dP)
    pb = pb - f / df
    return (pb + drho_g * H, pb)
end

# Newton for the bottom capillary pressure pb of a truncated plume:
# SP(pb + dP) − SP(pb) = SP_col, with dP = Δρ g H (plain values); co2lab's
# initial guess, iterated to round-off (co2lab stops at a residual of 1e-8
# times the largest dP of all columns)
function _bottom_capillary_pressure(m, SP_col, dP, pb)
    for _ in 1:200
        f = _interp_table(m.p, m.SP, pb) + SP_col - _interp_table(m.p, m.SP, pb + dP)
        # SP_col is at most SP's own scale, so this is round-off
        abs(f) <= 1e-13 * dP && return pb
        df = _interp_slope(m.p, m.SP, pb) - _interp_slope(m.p, m.SP, pb + dP)
        pb -= f / df
    end
    @warn "no convergence for the bottom capillary pressure of a truncated plume"
    return pb
end

# The free saturation the vertical model's hysteresis treatment gives
function _fringe_free_saturation(m::PScaledTable, sG, sGmax, rw, rg)
    m.hysteresis == :endpoint_scaling && return _free_saturation(sG, sGmax, rw, rg)
    error("'fixed residual' hysteresis for fringe fluids is roadmap Stage 6 item 2")
end

# ------------------------- Vertical-model interface ---------------------------

# Capillary heights Pc/(Δρ g) of the free saturation and of the maximum
# saturation (co2lab's upscaledSat2height, capillary case)
function _column_heights(m::PScaledTable, i, sG, sGmax, H, drho, rw, rg)
    drho_g = drho * Jutul.gravity_constant
    pc, _ = _fringe_capillary_pressure(m, _fringe_free_saturation(m, sG, sGmax, rw, rg),
                                       H, drho)
    pc_max, _ = _fringe_capillary_pressure(m, sGmax, H, drho)
    return (pc / drho_g, pc_max / drho_g)
end

# co2lab's krG_ptable and krW_simple (the fringe's brine relperm, a different
# function from the sharp interface's krW_simple: linear in the brine not
# bound as residual by the free CO2). No accretion layer: co2lab does not
# combine the fringe with rugosity.
function _column_relperm(m::PScaledTable, i, sG, sGmax, h, h_max, H, drho, dh,
                         krw, krg, rw, rg)
    drho_g = drho * Jutul.gravity_constant
    pc = h * drho_g
    pb = _max_c(pc - drho_g * H, 0.0)
    kr = _max_c(_interp_table(m.p, m.krP, pc), 0.0)
    # truncated plume: the part below the column's bottom
    pb > 0 && (kr -= _interp_table(m.p, m.krP, pb))
    krG = kr / (H * drho_g)

    sw = 1 - sG
    sw_eff = sw - _free_saturation(sG, sGmax, rw, rg) / (1 - rw) * rw
    return (_max_c(sw_eff, 0.0), krG)
end

_fine_saturations(::PScaledTable, i, zb, h, h_max, rw, rg) =
    error("reconstruction for fringe fluids is roadmap Stage 6 item 1b")

# The inventory's split of a fringe column: the mobile part by the current
# saturation profile, the immobile part by the maximum one (researcher,
# 2026-10-01); co2lab_compatible: co2lab's split by height (co2lab-issues.org)
function _category_saturations(m::PScaledTable, S, Smax, h, h_max, col, drho, rw, rg;
                               co2lab_compatible::Bool = false)
    co2lab_compatible && return _height_split(S, Smax, h, h_max, col, rw, rg)
    S_immob = Smax .* rg ./ (1 - rw)
    S_mob = S .- S_immob
    g = Jutul.gravity_constant
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
    SPat(z) = _interp_table(m.p, m.SP, max(pc - drho_g * clamp(z, 0.0, H), 0.0))
    total = SPat(0.0) - SPat(H)
    total > 0 || return 0.0
    return max(SPat(a) - SPat(b), 0.0) / total
end

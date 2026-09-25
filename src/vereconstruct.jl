# Reconstruction of fine-scale (3D) states from VE states (roadmap Stage 2
# item 6), as co2lab's VEstates23D does for its sharp-interface fluid:
#
# - pressure (pVE_to_p3D.m): the column's caprock brine pressure extended
#   hydrostatically with the brine density at that pressure,
#   p + ρ_w(p)·g·(z_cell − z_top), z_cell the 3D cell centroid depth;
# - saturation (height2finescaleSat.m, sharp_interface_h2s): from the plume
#   heights (h, h_max). Cells above h hold 1 − rw of CO2, the cell cut by h
#   the filled fraction of that; with residual CO2 (rg > 0) the zone between
#   h and h_max holds rg. The fill of a column is MRST's fillDegree.m.
#
# - dissolved CO2 (VEstates23D, rs part): the column's rs fills the brine from
#   the caprock down — saturated (dis_max) down to the depth where the brine
#   volume above holds the column's fraction rs/dis_max of it, none below.
#
# Capillary-fringe fluids (cap_fringe_h2s) are Stage 5.

# ------------------------------ Filling a column ------------------------------

"""
    _fill_degree(h, zb) -> (n, t)

MRST's `fillDegree` for one column with cumulative cell-bottom depths below
the caprock `zb` (increasing): `n` cells are completely filled by a plume of
height `h` (`zb[k] ≤ h`), and cell `n + 1` is filled to fraction `t` (0 when
all cells are filled). Negative `h` counts as 0.
"""
function _fill_degree(h, zb)
    h = max(h, zero(h))
    n = count(<=(h), zb)
    n == length(zb) && return (n, zero(h))
    # the two forms are MRST's, kept for identical rounding
    n == 0 && return (0, h / zb[1])
    return (n, -(zb[n] - h) / ((zb[n+1] - h) - (zb[n] - h)))
end

"""
    _column_saturations(zb, h, h_max, rw, rg) -> (s, smax)

CO2 saturation and maximum CO2 saturation of each cell of one column
(cumulative cell-bottom depths `zb`), from the plume height `h` and the
maximum plume height `h_max` — co2lab's `sharp_interface_h2s`. Returns two
new vectors, top cell first.
"""
function _column_saturations(zb, h, h_max, rw, rg)
    n, t = _fill_degree(h, zb)
    nk = length(zb)
    s = [k <= n ? 1 - rw : (k == n + 1 ? t * (1 - rw) : 0.0) for k in 1:nk]
    smax = copy(s)
    rg > 0 || return (s, smax)

    n_r, t_r = _fill_degree(h_max, zb)
    # the residual zone's partial cell is h's partial cell: only the part below h
    n_r == n && (t_r = max(t_r - t, 0.0))

    for k in 1:nk
        if k <= n_r && k > n + 1
            s[k], smax[k] = rg, 1 - rw
        elseif k <= n_r && k == n + 1
            s[k], smax[k] = s[k] + (1 - t) * rg, 1 - rw
        elseif k == n_r + 1
            s[k], smax[k] = s[k] + t_r * rg, s[k] + t_r * (1 - rw)
        end
    end
    return (s, smax)
end

"""
    _column_dissolved(zb, dz, brine, frac) -> Vector

Fraction of `dis_max` dissolved in each cell of one column (cumulative
cell-bottom depths `zb`, thicknesses `dz`, brine volume per unit height
`brine` = porosity × brine saturation) when the column holds the fraction
`frac` of its saturated amount: co2lab's `VEstates23D` finds the depth `h` at
which the brine above holds `frac` of the column's brine (Newton on
`integrateVertically`) and saturates everything above it. The integral is
piecewise linear, so `h` is found exactly here; with no brine at all, `h` is
co2lab's initial guess `frac·H`. Returns a new vector, top cell first.
"""
function _column_dissolved(zb, dz, brine, frac)
    target = frac * sum(brine .* dz)
    h = frac * zb[end]
    top, acc = 0.0, 0.0
    for k in eachindex(zb)
        v = brine[k] * dz[k]
        if v > 0 && acc + v >= target
            h = top + (target - acc) / brine[k]
            break
        end
        acc += v
        top = zb[k]
    end
    return first(_column_saturations(zb, h, h, 0.0, 0.0))
end

# ------------------------------- Reconstruction -------------------------------

"""
    ve_reconstruct(state, ve, columns, fluid; poro3D = nothing) -> NamedTuple

Fine-scale state on the 3D grid from a VE reservoir state, as co2lab's
`VEstates23D` for its sharp-interface fluid. `state` needs `:Pressure`,
`:sG`, `:sGmax` and, with dissolution, `:Rs` (a reservoir output state, or
the `:Reservoir` entry of a full state); `ve` is [`ve_domain`](@ref)'s
NamedTuple, `columns` the columns it was built on, `fluid` the
[`co2lab_ve_fluid`](@ref) of the run. `poro3D` (per 3D cell) weights the brine
volume for the dissolved CO2; `nothing` means uniform, as co2lab's default.

Returns the NamedTuple `(pressure, sG, sGmax, rs)` of new `Vector{Float64}`s,
one value per cell of the 3D mesh, in the mesh's cell order: brine pressure
extended hydrostatically from the caprock (Pa), CO2 saturation, maximum CO2
saturation, and dissolved CO2 ratio (`rs`; `nothing` without dissolution).
The plume heights are recomputed from `sG` and `sGmax` (co2lab's
`upscaledSat2height`). Arguments are only read.
"""
function ve_reconstruct(state, ve, columns::VEColumns, fluid::CO2LabVEFluid;
                        poro3D = nothing)
    geo, rw, rg = ve.geometry, fluid.res_water, fluid.res_gas
    p, sG, sGmax = state[:Pressure], state[:sG], state[:sGmax]
    nc = length(columns.cells)
    out = (pressure = zeros(nc), sG = zeros(nc), sGmax = zeros(nc))

    for col in 1:number_of_columns(columns)
        r = columns.pos[col]:(columns.pos[col+1] - 1)
        cells = columns.cells[r]

        h, h_max = _plume_heights(sG[col], sGmax[col], geo.H[col], rw, rg)
        s, smax = _column_saturations(view(geo.z, r), h, h_max, rw, rg)
        out.sG[cells], out.sGmax[cells] = s, smax

        # brine density at the caprock pressure, as co2lab's fluid.rhoW(pVE)
        rho_w = fluid.rho_ref[1] * _linear(1.0, fluid.rho_pvt[1], p[col])
        g = Jutul.gravity_constant
        out.pressure[cells] = p[col] .+ rho_w * g .* (geo.z_cell[r] .- geo.z_top[col])
    end
    fluid.dissolution && haskey(state, :Rs) || return (; out..., rs = nothing)

    # dissolved CO2 from the top of the column's brine down
    poro = isnothing(poro3D) ? ones(nc) : poro3D[columns.cells]
    rs = zeros(nc)
    for col in 1:number_of_columns(columns)
        r = columns.pos[col]:(columns.pos[col+1] - 1)
        brine = poro[r] .* (1 .- out.sG[columns.cells[r]])
        f = _column_dissolved(view(geo.z, r), view(geo.dz, r), brine,
                              state[:Rs][col] / fluid.dis_max)
        rs[columns.cells[r]] = f .* fluid.dis_max
    end
    return (; out..., rs = rs)
end

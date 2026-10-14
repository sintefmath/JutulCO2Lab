# Reconstruction of fine-scale (3D) states from VE states, as co2lab's
# VEstates23D:
#
# - pressure (co2lab's pVE_to_p3D): the column's caprock brine pressure
#   extended hydrostatically with the brine density at that pressure,
#   p + ρ_w(p)·g·(z_cell − z_top), z_cell the 3D cell centroid depth;
# - saturation: from the plume heights (h, h_max) by the fluid's vertical
#   model (`_fine_saturations`): co2lab's sharp_interface_h2s, or
#   cap_fringe_h2s for a capillary fringe. Sharp interface: cells above h hold
#   1 − rw of CO2, the cell cut by h its filled fraction of that, and with
#   residual CO2 (rg > 0) the zone between h and h_max holds rg;
# - dissolved CO2: the column's rs fills the brine from the caprock down,
#   saturated (dis_max) down to where the brine above holds the column's
#   fraction rs/dis_max of it, none below.

# ------------------------------ Filling a column ------------------------------

"""
    _fill_degree(h, zb) -> (n, t)

MRST's `fillDegree` for one column with cumulative cell-bottom depths below
the caprock `zb` (increasing). Returns `n`, the number of cells completely
filled by a plume of height `h` (`zb[k] ≤ h`), and `t`, the filled fraction
of cell `n + 1` (0 when all cells are filled). Negative `h` counts as 0.
"""
function _fill_degree(h, zb)
    h = max(h, zero(h))
    # n cells lie wholly above depth h (bottoms zb); t: the filled part of the next
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
        # wholly in the residual zone
        if k <= n_r && k > n + 1
            s[k], smax[k] = rg, 1 - rw
        # plume above, residual CO2 in the rest of the cell
        elseif k <= n_r && k == n + 1
            s[k], smax[k] = s[k] + (1 - t) * rg, 1 - rw
        # h_max's partial cell
        elseif k == n_r + 1
            s[k], smax[k] = s[k] + t_r * rg, s[k] + t_r * (1 - rw)
        end
    end
    return (s, smax)
end

"""
    _column_dissolved(zb, dz, brine, frac) -> Vector

Fraction of `dis_max` dissolved in each cell of one column holding the
fraction `frac` of its saturated amount, as co2lab's `VEstates23D`: the
brine is saturated from the caprock down to the depth `h` where the brine
above holds `frac` of the column's brine, and holds none below.

The column is given by cumulative cell-bottom depths `zb`, thicknesses `dz`
and brine volume per unit height `brine` (porosity × brine saturation).
`h` is found exactly (co2lab: Newton); with no brine at all it is co2lab's
initial guess `frac·H`. Returns a new vector, top cell first.
"""
function _column_dissolved(zb, dz, brine, frac)
    # brine to saturate: the fraction frac of the column's, filled from the top
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

# Modifies: `out.pressure`, `out.sG`, `out.sGmax` at `cell` (block `col`'s)
function _set_fine!(out, cell, p, sG, sGmax, col)
    out.pressure[cell], out.sG[cell], out.sGmax[cell] = p[col], sG[col], sGmax[col]
    return out
end

"""
    ve_reconstruct(state, ve, columns, fluid; poro3D = nothing) -> NamedTuple

The fine-scale state on the 3D grid from a VE reservoir state, as co2lab's
`VEstates23D`: for plotting and comparing VE results cell by cell.

# Arguments

- `state`: a VE reservoir state with `:Pressure`, `:sG`, `:sGmax` and, with
  dissolution, `:Rs` (a reservoir output state, or the `:Reservoir` entry of
  a full state).
- `ve`: [`ve_domain`](@ref)'s NamedTuple.
- `columns`: the columns `ve` was built on.
- `fluid`: the run's [`co2lab_ve_fluid`](@ref).
- `poro3D`: porosity per 3D cell, weighting the brine volume for the
  dissolved CO2; `nothing`: uniform, co2lab's default.

# Returns

The NamedTuple `(pressure, sG, sGmax, rs)` of new `Vector{Float64}`s, one
value per cell of the 3D mesh, in the mesh's cell order:

- `pressure`: brine pressure (Pa), extended hydrostatically from the caprock;
- `sG`, `sGmax`: CO2 saturation and maximum CO2 saturation;
- `rs`: dissolved CO2 ratio; `nothing` without dissolution.

# Notes

- The plume heights are recomputed from `sG` and `sGmax` by the fluid's
  vertical model (sharp interface: co2lab's `upscaledSat2height`).
- A fine 3D cell of a hybrid grid takes its block's own values.
"""
function ve_reconstruct(state, ve, columns::VEColumns, fluid::CO2LabVEFluid;
                        poro3D = nothing)
    geo, rw, rg = ve.geometry, fluid.res_water, fluid.res_gas
    p, sG, sGmax = state[:Pressure], state[:sG], state[:sGmax]
    rs_ve = haskey(state, :Rs) ? state[:Rs] : zeros(length(p))
    m = _vertical_model(fluid.vertical, ve.domain, fluid)
    # the columns' temperatures (tabulated properties), 30 °C where none
    d = ve.domain
    T = haskey(d, :temperature, Cells()) ? d[:temperature] :
        fill(Jutul.default_value(nothing, JutulDarcy.Temperature()), length(p))
    nc = length(columns.cells)
    out = (pressure = zeros(nc), sG = zeros(nc), sGmax = zeros(nc))

    for col in 1:number_of_columns(columns)
        r = columns.pos[col]:(columns.pos[col+1] - 1)
        cells = columns.cells[r]
        # a fine 3D cell of a hybrid grid: its own state (pressure at its
        # centroid)
        columns.is_ve[col] || (_set_fine!(out, only(cells), p, sG, sGmax, col); continue)

        drho = _density_difference(fluid, p[col], rs_ve[col], T[col])
        h, h_max = _column_heights(m, col, sG[col], sGmax[col], geo.H[col], drho, rw, rg)
        s, smax = _fine_saturations(m, col, view(geo.z, r), h, h_max, drho, rw, rg)
        out.sG[cells], out.sGmax[cells] = s, smax

        # brine density at the caprock pressure, as co2lab's fluid.rhoW(pVE)
        rho_w = fluid.rho_ref[1] * first(_shrinkage(fluid, p[col], T[col]))
        g = Jutul.gravity_constant
        out.pressure[cells] = p[col] .+ rho_w * g .* (geo.z_cell[r] .- geo.z_top[col])
    end
    fluid.dissolution && haskey(state, :Rs) || return (; out..., rs = nothing)

    # dissolved CO2 from the top of the column's brine down
    poro = isnothing(poro3D) ? ones(nc) : poro3D[columns.cells]
    rs = zeros(nc)
    for col in 1:number_of_columns(columns)
        r = columns.pos[col]:(columns.pos[col+1] - 1)
        columns.is_ve[col] || (rs[only(columns.cells[r])] = state[:Rs][col]; continue)
        brine = poro[r] .* (1 .- out.sG[columns.cells[r]])
        f = _column_dissolved(view(geo.z, r), view(geo.dz, r), brine,
                              state[:Rs][col] / fluid.dis_max)
        rs[columns.cells[r]] = f .* fluid.dis_max
    end
    return (; out..., rs = rs)
end

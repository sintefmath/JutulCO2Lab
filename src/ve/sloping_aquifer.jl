# co2lab's sloping aquifer (co2lab-ve/examples/helpers/sloping_aquifer.m with
# its test grid, make_testgrid.m): a tilted, slightly funnel-shaped aquifer
# with two Gaussian traps under its caprock, CO2 injected near the downdip
# end, the updip end open. Built here directly, without MRST.

# ---------------------------------- The grid ----------------------------------

"""
    sloping_aquifer_grid(; nx = 100, ny = 1, nz = 5, length = 2000.0,
                         width = 600.0, thickness = 20.0, bend = 5.0,
                         slope = 0.04, depth = 1000.0,
                         traps = ((0.5, 0.3, 20.0), (0.75, 0.2, 10.0)))
        -> Dict{String, Any}

co2lab's `make_testgrid` (without its noise and ripples) as a GeoEnergyIO
grid section: a Cartesian grid of `nx × ny × nz` cells, `length × width ×
thickness` m, whose node depths are shifted by the caprock surface `depth −
slope·x − bend·(y/W)(1 − y/W) − Σ traps`. Each trap `(xpos, w, height)` is a
Gaussian bump `height·exp(−d²/(w·W)²)`, centred at `x = xpos·L` on the
middle line, `d` the distance to it (along x only when `ny == 1`, co2lab's
cross-section). Node positions are computed as MRST's `linspace`, so they
match co2lab's grid exactly.

Returns a new `Dict` with `"cartDims"`, `"COORD"`, `"ZCORN"` and `"ACTNUM"`
(all active).
"""
function sloping_aquifer_grid(; nx::Integer = 100, ny::Integer = 1, nz::Integer = 5,
                              length = 2000.0, width = 600.0, thickness = 20.0,
                              bend = 5.0, slope = 0.04, depth = 1000.0,
                              traps = ((0.5, 0.3, 20.0), (0.75, 0.2, 10.0)))
    L, W = Float64(length), Float64(width)
    x, y = _mrst_linspace(0.0, L, nx + 1), _mrst_linspace(0.0, W, ny + 1)
    zl = _mrst_linspace(0.0, Float64(thickness), nz + 1)
    # the caprock depth at each pillar, in make_testgrid's order of operations
    function top(i, j)
        z = (depth - x[i] * slope) - bend * (y[j] / W) * (1 - y[j] / W)
        for (xpos, w, h) in traps
            d2 = (x[i] - xpos * L)^2 + (ny == 1 ? 0.0 : (y[j] - 0.5W)^2)
            z -= h * exp(-d2 / (w * W)^2)
        end
        return z
    end
    zt = [top(i, j) for i in 1:nx+1, j in 1:ny+1]

    # vertical pillars from the caprock to the bottom (x fastest, then y)
    coord = Float64[]
    for j in 1:ny+1, i in 1:nx+1
        append!(coord, (x[i], y[j], zt[i, j], x[i], y[j], zt[i, j] + zl[end]))
    end
    # ZCORN: per layer its top and bottom, each cell's corners on its pillars
    pillar(c) = c ÷ 2 + 1
    zcorn = [zl[k + t] + zt[pillar(ci), pillar(cj)]
             for ci in 1:2nx, cj in 1:2ny, t in 0:1, k in 1:nz]
    return Dict{String, Any}("cartDims" => (nx, ny, nz), "COORD" => coord,
                             "ZCORN" => vec(zcorn), "ACTNUM" => fill(true, nx, ny, nz))
end

# MRST's linspace(a, b, n): a + k·(b − a)/(n − 1), the last point b itself
_mrst_linspace(a, b, n) = [k == n - 1 ? b : a + k * ((b - a) / (n - 1)) for k in 0:n-1]

# --------------------------------- The case -----------------------------------

"""
    sloping_aquifer(; nz = 5, ny = 1, co2_rho = nothing, co2_mu = nothing,
                    kwargs...) -> NamedTuple

co2lab's sloping-aquifer case (`sloping_aquifer.m` with its defaults; `ny =
1` its cross-section, `'cross_sectional', true`), set up for JutulCO2Lab:
the grid of [`sloping_aquifer_grid`](@ref) (`kwargs` passed on) meshed by
GeoEnergyIO, its VE columns and [`ve_domain`](@ref) with co2lab's rock
(400 mD, porosity 0.2), and the case's data:

- brine 1050 kg/m³ and 0.8 mPa·s, the temperature 70 °C, the reference
  pressure `p_ref` (the mean hydrostatic pressure `1050·g·z` over the 3D
  cells) and rock compressibility 4.35e-5 /bar;
- CO2's density and viscosity at `p_ref` and 70 °C: by default from
  JutulDarcy's tables ([`jutuldarcy_property_tables`](@ref)); co2lab takes
  them from its own (`CO2props`), which `co2_rho` and `co2_mu` can give;
- the injector in pillar `(⌈nx/8⌉, ⌈2ny/5⌉)` (radius 0.1 m), 20 kt CO2 a
  year (`mass_rate`, kg/s);
- the open boundary: the updip end `x = L` (cross-section), or every
  lateral side; at hydrostatic brine pressure `1050·g·z`;
- the schedule: injection for 2 years in 24 steps (1 year in 12 for `ny >
  1`) after co2lab's ramp-up (`rampupTimesteps`, 7 halvings), then 2 years of
  migration in 12 steps.

Returns a new NamedTuple `(mesh, columns, ve, grid, sides, open_sides,
well_column, well_radius, rho_w, mu_w, co2_rho, co2_mu, temperature, p_ref,
rock_compressibility, mass_rate, dt, injecting)`: the mesh, `VEColumns`,
[`ve_domain`](@ref)'s result, the grid section, the boundary sides
([`ve_boundary_sides`](@ref)) and the indices of the open ones, the well's
column and radius (m), the constants (SI units; `rock_compressibility` 1/Pa
about `p_ref`), the report steps (s) and whether each injects.
"""
function sloping_aquifer(; nz::Integer = 5, ny::Integer = 1, co2_rho = nothing,
                         co2_mu = nothing, kwargs...)
    grid = sloping_aquifer_grid(; nz, ny, kwargs...)
    nx = grid["cartDims"][1]
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    columns, _ = ve_columns(mesh)
    nc = number_of_cells(mesh)
    ve = ve_domain(mesh, columns; perm_x = fill(400 * 1e-3 * _mrst_darcy(), nc),
                   poro = fill(0.2, nc))

    # co2lab's constants; its reference pressure over the 3D cells
    rho_w, g, T = 1050.0, Jutul.gravity_constant, 70 + 273.15
    z = tpfv_geometry(mesh).cell_centroids[3, :]
    p_ref = sum(rho_w * g .* z) / nc
    # CO2 at p_ref and T: JutulDarcy's unless given
    if isnothing(co2_rho) || isnothing(co2_mu)
        tab = jutuldarcy_property_tables()
        co2_rho = something(co2_rho, tab.density[2](p_ref, T))
        co2_mu = something(co2_mu, tab.viscosity[2](p_ref, T))
    end

    # the updip end open (cross-section), else every lateral side
    sides, _ = ve_boundary_sides(mesh, columns, ve.domain)
    dirs = ny == 1 ? (:right,) : (:left, :right, :lower, :upper)
    open_sides = ve_select_sides(sides, columns; directions = dirs)

    well = findfirst(==((cld(nx, 8), cld(2ny, 5))), columns.column_ij)
    # co2lab's schedule: its cross-section injects twice as long
    year = si_unit(:year)
    inj, n_inj = ny == 1 ? (2year, 24) : (year, 12)
    dt_inj = _rampup_timesteps(inj, inj / n_inj, 7)
    dt = [dt_inj; fill(2year / 12, 12)]
    injecting = [k <= length(dt_inj) for k in eachindex(dt)]
    return (mesh = mesh, columns = columns, ve = ve, grid = grid, sides = sides,
            open_sides = open_sides, well_column = well, well_radius = 0.1,
            rho_w = rho_w, mu_w = 8e-4, co2_rho = co2_rho, co2_mu = co2_mu,
            temperature = T, p_ref = p_ref, rock_compressibility = 4.35e-5 / 1e5,
            mass_rate = 20e6 / year, dt = dt, injecting = injecting)
end

# MRST's rampupTimesteps(time, dt, n): n + 1 steps doubling from dt/2ⁿ, then
# steps of dt, a final remainder merged into the last if tiny
function _rampup_timesteps(time, dt, n)
    time == 0 && return Float64[]
    init = [dt / 2^k for k in [n; n:-1:1]]
    any(>(time), cumsum(init)) && (init = init[cumsum(init) .< time])
    rem = fill(dt, floor(Int, (time - sum(init)) / dt))
    last_rem = isempty(rem) ? 0.0 : rem[end]
    final = time - sum(init) - sum(rem; init = 0.0)
    if final / last_rem <= 1e-6
        rem[end] = last_rem + final
        return [init; rem]
    end
    return [init; rem; final]
end

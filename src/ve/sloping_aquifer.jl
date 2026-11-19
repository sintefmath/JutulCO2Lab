# co2lab's sloping aquifer (its `sloping_aquifer` example and `make_testgrid`
# test grid): a tilted, slightly funnel-shaped aquifer with two Gaussian traps
# under its caprock, CO2 injected near the downdip end, the updip end open.
# Built here directly, without MRST.

# ---------------------------------- The grid ----------------------------------

"""
    sloping_aquifer_grid(; nx = 100, ny = 1, nz = 5, length = 2000.0,
                         width = 600.0, thickness = 20.0, bend = 5.0,
                         slope = 0.04, depth = 1000.0,
                         traps = ((0.5, 0.3, 20.0), (0.75, 0.2, 10.0)))
        -> Dict{String, Any}

The sloping aquifer's grid: co2lab's `make_testgrid` (without its noise and
ripples) as a GeoEnergyIO grid section. Node positions follow MRST's
`linspace`, so they match co2lab's grid exactly.

# Arguments

- `nx`, `ny`, `nz`: cells along x, y and down; `ny = 1` is co2lab's
  cross-section.
- `length`, `width`, `thickness`: the box's size `L × W × H` (m).
- `depth`, `slope`, `bend`: the caprock surface before the traps,
  `depth − slope·x − bend·(y/W)(1 − y/W)` (m).
- `traps`: tuples `(xpos, w, height)`, each subtracting a Gaussian bump
  `height·exp(−d²/(w·W)²)` (m) centred at `x = xpos·L` on the middle line,
  `d` the distance to it (along x only when `ny == 1`).

# Returns

A new `Dict` with `"cartDims"`, `"COORD"`, `"ZCORN"` and `"ACTNUM"` (all
active).
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

co2lab's sloping-aquifer case with its defaults, set up for JutulCO2Lab: the
mesh, its VE columns and domain, and the case's data, ready for a model and a
schedule.

# Arguments

- `nz`, `ny`: cells down and across; `ny = 1` is co2lab's cross-section
  (`'cross_sectional', true`).
- `co2_rho`, `co2_mu`: CO2's density (kg/m³) and viscosity (Pa·s) at `p_ref`
  and 70 °C; default from JutulDarcy's tables
  ([`jutuldarcy_property_tables`](@ref)). co2lab takes them from its own
  `CO2props`; pass those for comparisons.
- `kwargs...`: passed to [`sloping_aquifer_grid`](@ref).

# Returns

A new NamedTuple (SI units):

- `mesh`, `columns`, `ve`, `grid`: the GeoEnergyIO mesh, its `VEColumns`,
  [`ve_domain`](@ref)'s result (co2lab's rock: 400 mD, porosity 0.2) and the
  grid section;
- `sides`, `open_sides`: the boundary sides ([`ve_boundary_sides`](@ref)) and
  the indices of the open ones: the updip end `x = L` for the cross-section,
  every lateral side otherwise; held at hydrostatic pressure `1050·g·z`;
- `well_column`, `well_radius`: the injector's column, in pillar
  `(⌈nx/8⌉, ⌈2ny/5⌉)`, and its radius (0.1 m);
- `rho_w`, `mu_w`: brine density and viscosity (1050 kg/m³, 0.8 mPa·s);
- `co2_rho`, `co2_mu`: as above;
- `temperature`: 70 °C, in K;
- `p_ref`: the reference pressure, the mean hydrostatic pressure `1050·g·z`
  over the 3D cells;
- `rock_compressibility`: 4.35e-5 /bar, in 1/Pa, about `p_ref`;
- `mass_rate`: the injection rate, 20 kt CO2 a year, in kg/s;
- `dt`: the report steps (s): injection for 2 years in 24 steps (1 year in
  12 for `ny > 1`) after co2lab's ramp-up ([`rampup_timesteps`](@ref), 7
  halvings), then 2 years of migration in 12 steps;
- `injecting`: whether each report step injects.
"""
function sloping_aquifer(; nz::Integer = 5, ny::Integer = 1, co2_rho = nothing,
                         co2_mu = nothing, kwargs...)
    grid = sloping_aquifer_grid(; nz, ny, kwargs...)
    nx = grid["cartDims"][1]
    mesh = _mesh_from_grid_section(grid)
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
    dt_inj = rampup_timesteps(inj, inj / n_inj, 7)
    dt = [dt_inj; fill(2year / 12, 12)]
    injecting = [k <= length(dt_inj) for k in eachindex(dt)]
    return (mesh = mesh, columns = columns, ve = ve, grid = grid, sides = sides,
            open_sides = open_sides, well_column = well, well_radius = 0.1,
            rho_w = rho_w, mu_w = 8e-4, co2_rho = co2_rho, co2_mu = co2_mu,
            temperature = T, p_ref = p_ref, rock_compressibility = 4.35e-5 / 1e5,
            mass_rate = 20e6 / year, dt = dt, injecting = injecting)
end

"""
    rampup_timesteps(time, dt, n) -> Vector{Float64}

MRST's `rampupTimesteps`: time steps covering `time` that start small and
grow to `dt`, for a gentle start of an injection.

# Arguments

- `time`: the total time.
- `dt`: the target step length (same unit as `time`).
- `n`: the number of halvings in the ramp-up.

# Returns

A new vector of step lengths (the unit of `time`); empty for `time = 0`:

- `n + 1` steps doubling from `dt/2ⁿ` (cut short if `time` is reached
  first);
- then steps of `dt`;
- then a shorter step for the remainder, merged into the last step if below
  `1e-6 dt`.
"""
function rampup_timesteps(time, dt, n)
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

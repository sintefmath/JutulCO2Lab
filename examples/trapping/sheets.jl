# Trap analysis per sheet: an aquifer with an internal sealing layer. CO2
# injected below the layer pools in the layer's own trap, spills over the
# layer's rim and rises to the caprock above, where it fills the caprock's
# traps. The grid's VE columns form two sheets — the caprock, and the surface
# under the layer — each with its own trap analysis (ve_sheets,
# trap_surface(...; sheet), sheet_spill_targets).
#
#   julia --project=examples examples/trapping/sheets.jl
#
# Self-contained (the grid is built here). CO2 and brine properties are
# JutulDarcy's. Plotting backend: GLMakie (a window per figure; 3D views
# rotatable with the mouse; close it to continue) by default;
# JUTULCO2LAB_BACKEND=CairoMakie writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, JutulDarcy, Printf

year = si_unit(:year)

# ------------------------------- 1. The aquifer -------------------------------

# co2lab's test aquifer (make_testgrid), 4 x 2 km, 40 m thick, 40 x 20 x 20
# cells: the caprock rises 1 % towards x = 4 km, with two domes (traps) at
# x = 2 and 3 km.
nx, ny, nz = 40, 20, 20
grid = sloping_aquifer_grid(; nx, ny, nz, length = 4000.0, width = 2000.0,
                            thickness = 40.0, slope = 0.01,
                            traps = ((0.5, 0.3, 20.0), (0.75, 0.2, 10.0)))
mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]

# The sealing layer: the faces between layers 10 and 11 under a 1.3 x 1 km
# patch around the first dome (x 1.3-2.6 km, y 0.5-1.5 km). It follows the
# caprock's shape, so it has a dome (a trap) of its own.
patch(i, j) = 14 <= i <= 26 && 6 <= j <= 15
N = get_neighborship(mesh)
function in_layer(f)
    a, b = ijk[N[1, f]], ijk[N[2, f]]
    return a[1:2] == b[1:2] && patch(a[1], a[2]) && minmax(a[3], b[3]) == (10, 11)
end
layer = filter(in_layer, axes(N, 2))
# the layer seals: transmissibility multiplier 0 on its faces
mult = ones(size(N, 2))
mult[layer] .= 0.0

# VE columns, cut at the layer: below it, a second VE column in each pillar
columns, _ = ve_columns(mesh; transmissibility_multiplier = mult, below_holes = :stack)
ve = ve_domain(mesh, columns; perm_x = fill(500e-15, nz * nx * ny),
               poro = fill(0.2, nz * nx * ny))

# ---------------------------- 2. Sheets and traps -----------------------------

# Two sheets: the caprock (the first VE column of every pillar) and the
# surface under the layer. Each gets its own trap analysis.
sheets = ve_sheets(mesh, columns)
@printf("%d VE columns in %d sheets: %s\n", number_of_columns(columns),
        sheets.nsheets, join(["$(count(==(k), sheets.sheet))" for k in 1:2], " + "))
surfaces = [trap_surface(mesh, columns; sheet = k) for k in 1:sheets.nsheets]
results = [trap_analysis(s) for s in surfaces]
for (k, (s, r)) in enumerate(zip(surfaces, results))
    v = trap_volumes(s, r; poro = fill(0.2, length(s.blocks)))
    @printf("sheet %d: %d traps, pore volume %s m³\n", k, length(v),
            join([@sprintf("%.0f", x) for x in v], ", "))
end

# Where CO2 goes that leaves the lower sheet: over its rim, into the
# neighbouring pillar and up to the column it reaches there (on the caprock)
spills = sheet_spill_targets(mesh, columns, surfaces[2])
@printf("lower sheet: %d rim edges, all leading to sheet %s\n", length(spills),
        join(unique(x.target_sheet for x in spills), ", "))

# The rim edges the lower trap's spill paths (rivers) end at
river_ends = Set(last(r) for rs in results[2].rivers for r in rs)
spill_out = filter(x -> x.polygon in river_ends, spills)

# ------------------------------ 3. The traps in 3D -----------------------------

# The two sheets side by side, at the same scale: the caprock with the points
# where CO2 spilled from the layer arrives (red), and the surface under the
# layer with the rim edges its trap spills over (red)
s1, s2 = surfaces
function sheet_axis(pos, title)
    ax = ve_axis3(pos, s1; exaggeration = 8, view = (-50, 35), title = title)
    # both panels at the caprock's extent; depth bounds deepest first, as
    # limits! with them shallowest first would turn the depth axis upwards
    limits!(ax, extrema(s1.x)..., extrema(s1.y)..., reverse(extrema([s1.z; s2.z]))...)
    return ax
end
poly_of = Dict(b => i for (i, b) in enumerate(s1.blocks))
arrive = unique(poly_of[x.target] for x in spill_out)
rim = [Point3f(s2.x[n], s2.y[n], s2.z[n] - 0.5) for x in spill_out
       for n in s2.edges[:, x.edge]]

fig = Figure(size = (1600, 700))
ax = sheet_axis(fig[1, 1], "Caprock: traps; red: CO2 arriving from the layer")
plot_traps!(ax, s1, results[1])
scatter!(ax, Point3f.(s1.cx[arrive], s1.cy[arrive], s1.zc[arrive] .- 0.5);
         color = :red, markersize = 14)
ax = sheet_axis(fig[1, 2], "Under the layer: its trap; red: where it spills over")
plot_traps!(ax, s2, results[2])
linesegments!(ax, rim; color = :red, linewidth = 4)
show_figure(fig, "sheets_traps")

# ---------------------------- 4. Inject and migrate ----------------------------

# 10 kt of CO2 a year for 20 years below the layer, upslope of the patch's
# downdip edge; then 280 years of migration. JutulDarcy's CO2-brine tables at
# 40 °C; residual saturations 0.2 (hysteresis). The updip end (x = 4 km) is
# open, the other sides closed.
fluid = co2lab_ve_fluid(; residual = (0.2, 0.2))
injector = only(b for b in 1:number_of_columns(columns)
                if columns.column_ij[b] == (17, 10) && sheets.sheet[b] == 2)
sides, _ = ve_boundary_sides(mesh, columns, ve.domain)
open_end = [d == 2 && columns.column_ij[c][1] == nx ? :open : :closed
            for (c, d) in zip(sides.column, sides.direction)]
case = setup_injection_case((mesh = mesh, columns = columns, ve = ve),
                            [(injector, 10e6 / year)]; fluid,
                            boundary = open_end, temperature = 273.15 + 40,
                            inj_time = 20year, inj_steps = 10,
                            mig_time = 280year, mig_steps = 28)
t = @elapsed result = simulate_injection_case(case; info_level = -1)
states = result.states
@printf("300 years simulated in %.0f s\n", t)

# CO2 in each sheet over time (kg): the column masses of each sheet
mass(st, k) = sum(st[:TotalMasses][2, sheets.sheet .== k])
for k in (10, 15, 20, 38)
    @printf("year %3.0f: %5.1f kt under the layer, %5.1f kt under the caprock\n",
            sum(case.dt[1:k]) / year, mass(states[k], 2) / 1e6, mass(states[k], 1) / 1e6)
end

# ------------------------------ 5. The plume in 3D -----------------------------

# The plume height on each sheet, side by side on one colour scale: at the
# end of injection, and at 300 years
function plume_view(k, when)
    h = states[k][:PlumeHeights][1, :]
    hmax = max(maximum(h), 1e-3)
    # one colour scale for both sheets: the colours computed here
    colours(s) = [cgrad(:plasma)[v / hmax] for v in h[s.blocks]]
    fig = Figure(size = (1600, 700))
    for (n, (s, name)) in enumerate(((s1, "under the caprock"), (s2, "under the layer")))
        ax = sheet_axis(fig[1, n], "$when: CO2 $name")
        plot_ve_surface!(ax, s; color = colours(s), strokewidth = 0)
    end
    Colorbar(fig[1, 3]; limits = (0, hmax), colormap = :plasma,
             label = "CO2 plume height (m)")
    return fig
end
show_figure(plume_view(10, "20 years (end of injection)"), "sheets_plume_20")
show_figure(plume_view(length(states), "300 years"), "sheets_plume_300")

# --------------------------- 6. A cross-section -------------------------------

# The reconstructed CO2 saturation on the vertical section y = 1 km (j = 10)
# through the injector, at four times: CO2 pools under the layer (the dashed
# line), spills at its updip rim and rises to the caprock.
section = findall(x -> x[2] == 10, ijk)
function cell_nodes(c)
    f, b = mesh.faces, mesh.boundary_faces
    return Set([(n for g in f.cells_to_faces[c] for n in f.faces_to_nodes[g])...,
                (n for g in b.cells_to_faces[c] for n in b.faces_to_nodes[g])...])
end
function outline(c)
    pts = [Point2f(p[1], p[3]) for p in mesh.node_points[collect(cell_nodes(c))]]
    m = sum(pts) / length(pts)
    # the four corners seen from the side: the two leftmost and two rightmost
    left = sort(filter(p -> p[1] < m[1], pts); by = p -> p[2])
    right = sort(filter(p -> p[1] > m[1], pts); by = p -> p[2])
    return Makie.Polygon([left[1], right[1], right[end], left[end]])
end
shapes = outline.(section)
# the layer on the section: its faces' centroids, in x order
fc = tpfv_geometry(mesh).face_centroids
on_section = sort(filter(f -> ijk[N[1, f]][2] == 10, layer); by = f -> fc[1, f])
layer_line = [Point2f(fc[1, f], fc[3, f]) for f in on_section]
fig = Figure(size = (1200, 900))
for (n, k) in enumerate((4, 10, 16, 38))
    sG = ve_reconstruct(states[k], case.ve, columns, fluid).sG
    local ax = Axis(fig[n, 1]; title = @sprintf("%.0f years", sum(case.dt[1:k]) / year),
                    ylabel = "depth (m)", xlabel = n == 4 ? "x (m)" : "",
                    yreversed = true)
    poly!(ax, shapes; color = sG[section], colorrange = (0, 0.8), colormap = :plasma)
    lines!(ax, layer_line; color = :white, linestyle = :dash, linewidth = 2)
end
Colorbar(fig[:, 2]; colorrange = (0, 0.8), colormap = :plasma, label = "CO2 saturation")
show_figure(fig, "sheets_section")

# Every report step in Jutul's interactive 3D viewer (GLMakie only): the
# slider steps through time
if backend == :GLMakie
    sat = [Dict(:CO2_saturation => ve_reconstruct(s, case.ve, columns, fluid).sG)
           for s in states]
    show_figure(plot_interactive(mesh, sat; new_window = false), "sheets_viewer")
end

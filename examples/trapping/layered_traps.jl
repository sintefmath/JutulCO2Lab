# Traps under internal layers, after Møyner and Nilsen (2019), Fig. 7: a
# tilted cross-section with an undulating caprock and three partial sealing
# layers. CO2 injected at the lower left fills the traps under the lowest
# layer, spills over its up-dip end, rises to the next layer, and so on up to
# the caprock; the right end is open.
#
# Two separate parts:
# - the spill-point analysis, on the VE grid's sheets (one per layer and the
#   caprock): traps, spill paths, where each sheet's CO2 goes when it spills;
# - the simulation, VE and 3D side by side, from one fine-scale fluid, with
#   run times.
#
#   julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/trapping/layered_traps.jl
#
# Self-contained (the grid is built here). CO2 and brine properties are
# JutulDarcy's. The 3D cells are 1 m tall, so that the plumes (several metres
# under the layers) are resolved and the 3D run's relperm smearing is small.
# Plotting backend: GLMakie (a window per figure; close it to continue) by
# default; JUTULCO2LAB_BACKEND=CairoMakie writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, JutulDarcy, Printf

year = si_unit(:year)

# ------------------------------- 1. The aquifer -------------------------------

# 1000 m long, 100 m thick, 100 m wide (one cell); 200 columns of 100 cells.
# The caprock rises 2.5° to the right, with undulations of three wavelengths;
# every layer below follows it.
nx, nz, L, H, W = 200, 100, 1000.0, 100.0, 100.0
top(x) = 1000.0 - tand(2.5) * x - 7sin(2π * x / 380) - 3sin(2π * x / 150 + 1) -
         1.5sin(2π * x / 60 + 2)
function aquifer_mesh()
    coord = Float64[]
    for j in 0:1, i in 0:nx
        x = i * L / nx
        append!(coord, [x, j * W, 0.0, x, j * W, 2000.0])
    end
    # each cell's corners on its pillars: depth top(x) + (k − 1) H/nz
    zc = [top((ii ÷ 2) * L / nx) + (k - 1 + s) * H / nz
          for k in 1:nz for s in (0, 1) for jj in 1:2 for ii in 1:2nx]
    grid = Dict("cartDims" => [nx, 1, nz], "COORD" => coord, "ZCORN" => zc,
                "ACTNUM" => fill(true, nx, 1, nz))
    return JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
end
mesh = aquifer_mesh()
nc = number_of_cells(mesh)
ijk = [cell_ijk(mesh, c) for c in 1:nc]
perm, poro = fill(100 * 1e-3 * si_unit(:darcy), nc), fill(0.2, nc)

# Three sealing layers (no flow across them), each between two cell layers
# over part of the length, staggered up-dip: deep on the left, shallow on
# the right. (layer below which they lie, columns)
layers = [(70, 1:95), (45, 55:150), (20, 115:190)]
N = get_neighborship(mesh) # cell neighbors for each internal face
function on_layer(f, (k, cols))
    a, b = ijk[N[1, f]], ijk[N[2, f]]
    return a[1] == b[1] && a[1] in cols && minmax(a[3], b[3]) == (k, k + 1)
end
sealed = [f for f in axes(N, 2) if any(l -> on_layer(f, l), layers)]

# --------------------------- 2. Spill-point analysis --------------------------

# The VE grid: the pillars cut at the layers; each piece a VE column.
columns, _ = ve_columns(mesh; sealed_faces = sealed, below_holes = :stack)
sheets = ve_sheets(mesh, columns)
# Closed as in the simulation: the section's sides y = 0 and y = W and its
# left end (sealed edges, which the trap analysis and the spill targets treat
# as closed); open are only the right end and the layers' rims.
function closed_sides(s)
    on(k, f, v) = all(n -> f(s, n) == v, s.edges[:, k])
    closed(k) = on(k, (s, n) -> s.y[n], 0.0) || on(k, (s, n) -> s.y[n], W) ||
                on(k, (s, n) -> s.x[n], 0.0)
    return (; s..., sealed = [s.sealed[k] || (s.edge_columns[2, k] == 0 && closed(k))
                              for k in axes(s.edges, 2)])
end
surfaces = [closed_sides(trap_surface(mesh, columns; sheet = k))
            for k in 1:sheets.nsheets]
results = [trap_analysis(s) for s in surfaces]
spills = [sheet_spill_targets(mesh, columns, s) for s in surfaces]

# a sheet by what lies above it: the caprock, or a layer (its x range)
function sheet_name(k)
    sheets.top[k] && return "caprock"
    a, b = extrema(first.(columns.column_ij[surfaces[k].blocks]))
    return @sprintf("layer at x %.0f-%.0f m", (a - 1) * L / nx, b * L / nx)
end

# Each sheet: its depth range, traps and their pore volumes, and where CO2
# spilling from it goes
for (k, (s, r)) in enumerate(zip(surfaces, results))
    v = trap_volumes(s, r; poro = fill(0.2, length(s.blocks)))
    to = unique(x.target_sheet for x in spills[k] if x.target > 0)
    @printf("sheet %d: columns %3d-%3d, %d traps, %6.0f m³ of pore volume; ",
            k, extrema(first.(columns.column_ij[s.blocks]))..., length(v), sum(v))
    println("spills to ", isempty(to) ? "the boundary" : "sheet " * join(to, ", "))
end

# ------------------- 3. Figure: the spill-point analysis ----------------------

# The section at y = 0: the aquifer's outline, the layers (black), each
# sheet's traps filled down to their spill depth (one colour per sheet), the
# spill paths along the sheets (rivers, thin), and where CO2 spilling from a
# layer's trap goes up (red: from the end of its spill path, over the rim, up
# to the sheet above)
xs = range(0, L; length = nx + 1)
# a column's two top corners at y = 0: its top edge on the section
function top_edge(s, p)
    n = filter(n -> s.y[n] == 0.0, unique(s.polygons[p]))
    a, b = sort(n; by = n -> s.x[n])
    return (Point2f(s.x[a], s.z[a]), Point2f(s.x[b], s.z[b]))
end
# a trap column filled from its top edge down to the trap's spill depth
function trap_fill(s, p, z)
    a, b = top_edge(s, p)
    return Makie.Polygon([a, b, Point2f(b[1], max(b[2], z)),
                          Point2f(a[1], max(a[2], z))])
end
centre(s, p) = Point2f(s.cx[p], s.zc[p])
colours = Makie.wong_colors()

fig = Figure(size = (1400, 600))
ax = Axis(fig[1, 1]; title = "Spill-point analysis per sheet: traps filled to " *
          "their spill depth; red: spill to the sheet above",
          xlabel = "x (m)", ylabel = "depth (m)", yreversed = true)
lines!(ax, xs, top.(xs); color = :gray40)
lines!(ax, xs, top.(xs) .+ H; color = :gray40)
for (k, cols) in layers
    x = xs[first(cols):last(cols)+1]
    lines!(ax, x, top.(x) .+ k * H / nz; color = :black, linewidth = 3)
end
for (k, (s, r)) in enumerate(zip(surfaces, results))
    c = colours[mod1(k, length(colours))]
    trapped = findall(>(0), r.traps)
    poly!(ax, [trap_fill(s, p, r.trap_z[r.traps[p]]) for p in trapped]; color = (c, 0.8))
    rivers = [centre(s, p) for rs in r.rivers for rv in rs
              for (a, b) in zip(rv[1:end-1], rv[2:end]) for p in (a, b)]
    isempty(rivers) || linesegments!(ax, rivers; color = c, linewidth = 1.5)
    # spills to another sheet from where the spill paths end (not the ends)
    ends = Set(last(rv) for rs in r.rivers for rv in rs)
    for x in spills[k]
        (x.polygon in ends && x.target_sheet > 0) || continue
        t = findfirst(==(x.target), surfaces[x.target_sheet].blocks)
        from, to = centre(s, x.polygon), centre(surfaces[x.target_sheet], t)
        lines!(ax, [from, Point2f(from[1], to[2])]; color = :red, linewidth = 2)
        scatter!(ax, [Point2f(from[1], to[2])]; color = :red, marker = :utriangle)
    end
end
scatter!(ax, [Point2f(2.5L / nx, top(0) + H - 5)]; color = :black, markersize = 14,
         marker = :star5)
text!(ax, Point2f(2.5L / nx + 10, top(0) + H - 5); text = "injector", fontsize = 12,
      align = (:left, :center))
Legend(fig[1, 2], [PolyElement(color = colours[mod1(k, length(colours))])
                   for k in 1:sheets.nsheets], sheet_name.(1:sheets.nsheets), "sheet")
show_figure(fig, "layered_traps_analysis")

# ---------------------------- 4. VE and 3D runs -------------------------------

# One fine-scale fluid for both: no capillary pressure, linear relperms (the
# VE columns then hold a sharp interface), no residual saturations (CO2 is
# trapped structurally only, so the run can be checked against the
# spill-point analysis); JutulDarcy's CO2–brine tables at 40 °C. 4 kt of CO2 a
# year for 30 years at the lower left (column 3, its bottom 10 m), then 970
# years of migration; the right end open at hydrostatic pressure.
fine_fluid = FineScaleFluid(; rock_types = [RockCurves(; kr3D = 1)])
dt = [fill(year, 30); fill(10year, 97)]
rate = 4e6 / year

# The run on the VE columns (`kind = :VE`) or on the 3D cells (`:3D`), for the
# report steps `dt`
function setup_run(kind, dt)
    cols = kind == :VE ? columns :
           first(ve_columns(mesh; sealed_faces = sealed, is_ve = fill(false, nc),
                            below_holes = :stack))
    ve = ve_domain(mesh, cols; perm_x = perm, poro, temperature = 273.15 + 40)
    fluid = co2lab_ve_fluid(fine_fluid, mesh, cols; poro, perm_x = perm)
    model = setup_co2lab_model(ve.domain; wells = [injector(kind, ve, cols)], fluid)
    p = ve_hydrostatic_profile(fluid, ve)
    sides, _ = ve_boundary_sides(mesh, cols, ve.domain)
    bc = ve_boundary_conditions(ve, sides, ve_select_sides(sides, cols;
                                directions = (:right,)); pressure = p)
    inject = co2_injector_control(TotalRateTarget(rate / fluid.rho_ref[2]), fluid)
    on = setup_reservoir_forces(model; control = Dict(:inj => inject), bc)
    off = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()), bc)
    forces = [t <= 30year ? on : off for t in cumsum(dt)]
    state0 = setup_co2lab_state(model; pressure = p.(ve.domain[:reference_depth]),
                                forces = on)
    case = JutulCase(model, dt, forces; state0, parameters = setup_parameters(model))
    return (case = case, ve = ve, columns = cols, fluid = fluid)
end

# The injector: the VE column at the bottom of column 3, or its lowest ten 3D
# cells (Peaceman's index, a small wellbore that starts full of CO2)
function injector(kind, ve, cols)
    cell(k) = findfirst(==((3, 1, k)), ijk)
    kind == :VE && return setup_ve_well(ve, cols.partition[cell(nz)]; name = :inj)
    perf = [cols.partition[cell(k)] for k in nz-9:nz]
    WI = [JutulDarcy.compute_peaceman_index(JutulDarcy.cell_dims(mesh, cell(k)),
                                            perm[cell(k)], 0.1) for k in nz-9:nz]
    return setup_well(ve.domain, perf; name = :inj, WI, simple_well = true,
                      volume_multiplier = 1e-3)
end

# Each run timed after a two-step run of the same kind (Julia compiles then)
function timed_run(kind)
    warm = setup_run(kind, dt[1:2])
    sim, cfg = setup_co2lab_simulator(warm.case; info_level = 1)
    simulate_reservoir(warm.case; simulator = sim, config = cfg)
    r = setup_run(kind, dt)
    sim, cfg = setup_co2lab_simulator(r.case; info_level = 1)
    t = @elapsed res = simulate_reservoir(r.case; simulator = sim, config = cfg)
    steps = sum(length(rep[:ministeps]) for rep in res.result.reports)
    return (; r..., states = res.states, seconds = t, steps = steps)
end
runs = [kind => timed_run(kind) for kind in (:VE, Symbol("3D"))]
for (kind, r) in runs
    @printf("%-3s %6d cells: 1000 years in %6.1f s (%d time steps)\n", kind,
            number_of_columns(r.columns), r.seconds, r.steps)
end

# ------------------------- 5. Where the CO2 ends up ---------------------------

# CO2 per sheet: the VE columns of each sheet, and the 3D cells in them
sheet_of_cell = [sheets.sheet[columns.partition[c]] for c in 1:nc]
function co2_per_sheet(r, st)
    m = zeros(sheets.nsheets)
    for b in 1:number_of_columns(r.columns), c in column_cells(r.columns, b)
        # a VE column's mass on its sheet; a 3D cell (one per block) on its own
        n = length(column_cells(r.columns, b))
        m[sheet_of_cell[c]] += st[:TotalMasses][2, b] / n
    end
    return m
end
times = cumsum(dt)
ve_run, fine_run = last(runs[1]), last(runs[2])
injected = rate .* min.(times, 30year)

# Per trap of the spill-point analysis: the 3D cells inside it (in its
# columns, above its spill depth) and the fraction of their pore volume that
# holds CO2 after 1000 years, VE and 3D: about 1 for a trap on the CO2's
# path, 0 for one it does not reach
depth = tpfv_geometry(mesh).cell_centroids[3, :]
pv = poro .* tpfv_geometry(mesh).volumes
function trap_cells(k, t)
    s, r = surfaces[k], results[k]
    return [c for p in findall(==(t), r.traps)
            for c in column_cells(columns, s.blocks[p]) if depth[c] < r.trap_z[t]]
end
final_sG = [ve_reconstruct(r.states[end], r.ve, r.columns, r.fluid).sG
            for (_, r) in runs]
filled(sG, cells) = sum(sG[cells] .* pv[cells]) / sum(pv[cells])
@printf("\n%-22s %5s %12s %10s %8s %8s\n", "sheet", "trap", "x (m)", "pore (m³)",
        "VE", "3D")
for k in 1:sheets.nsheets, t in eachindex(results[k].trap_z)
    cells = trap_cells(k, t)
    isempty(cells) && continue
    x = extrema(tpfv_geometry(mesh).cell_centroids[1, cells])
    @printf("%-22s %5d %5.0f-%5.0f %10.0f %8.2f %8.2f\n", sheet_name(k), t, x...,
            sum(pv[cells]), filled(final_sG[1], cells), filled(final_sG[2], cells))
end

# CO2 per sheet after 1000 years, and what has left through the open end
@printf("\n%-22s %8s %8s\n", "CO2 after 1000 years", "VE (kt)", "3D (kt)")
in_place = [co2_per_sheet(r, r.states[end]) ./ 1e6 for (_, r) in runs]
for k in 1:sheets.nsheets
    @printf("%-22s %8.1f %8.1f\n", sheet_name(k), in_place[1][k], in_place[2][k])
end
@printf("%-22s %8.1f %8.1f   (of %.0f kt injected)\n", "left the domain",
        injected[end] / 1e6 - sum(in_place[1]), injected[end] / 1e6 - sum(in_place[2]),
        injected[end] / 1e6)

# ----------------------- 6. Figures: VE and 3D results -------------------------

# The CO2 saturation on the section (VE columns reconstructed on the 3D
# cells), with the layers (black) and each trap's spill depth from the
# spill-point analysis (white): at the end of injection and after 1000 years
function cell_outline(c)
    bf = mesh.boundary_faces
    nodes = Set(n for f in bf.cells_to_faces[c] for n in bf.faces_to_nodes[f])
    pts = [Point2f(p[1], p[3]) for p in mesh.node_points[collect(nodes)] if p[2] == 0]
    m = sum(pts) / length(pts)
    return Makie.Polygon(sort(pts; by = p -> atan(p[2] - m[2], p[1] - m[1])))
end
outlines = cell_outline.(1:nc)
# per trap a line at its spill depth, across its columns
function spill_line(s, r, t)
    x = [q[1] for p in findall(==(t), r.traps) for q in top_edge(s, p)]
    return (Point2f(minimum(x), r.trap_z[t]), Point2f(maximum(x), r.trap_z[t]))
end
spill_lines = [q for (s, r) in zip(surfaces, results) for t in eachindex(r.trap_z)
               for q in spill_line(s, r, t)]

fig = Figure(size = (1600, 900))
for (row, (kind, r)) in enumerate(runs), (col, k) in enumerate((30, length(dt)))
    local ax = Axis(fig[row, col]; title = @sprintf("%s, %.0f years (%.0f s)", kind,
                                                    times[k] / year, r.seconds),
                    xlabel = "x (m)", ylabel = "depth (m)", yreversed = true)
    sG = ve_reconstruct(r.states[k], r.ve, r.columns, r.fluid).sG
    poly!(ax, outlines; color = sG, colorrange = (0, 1), colormap = :plasma)
    linesegments!(ax, spill_lines; color = :white, linewidth = 2)
    for (kk, cols) in layers
        x = xs[first(cols):last(cols)+1]
        lines!(ax, x, top.(x) .+ kk * H / nz; color = :black, linewidth = 2)
    end
end
Colorbar(fig[:, 3]; colorrange = (0, 1), colormap = :plasma, label = "CO2 saturation")
show_figure(fig, "layered_traps_saturation")

# CO2 under each sheet over time, VE (solid) and 3D (dashed), and what has
# left through the open end
fig = Figure(size = (1000, 600))
ax = Axis(fig[1, 1]; xlabel = "years", ylabel = "CO2 (kt)", xscale = log10,
          title = "CO2 per sheet: VE (solid) and 3D (dashed)")
for k in 1:sheets.nsheets
    c = colours[mod1(k, length(colours))]
    for (kind, r) in runs
        m = [co2_per_sheet(r, st)[k] for st in r.states] ./ 1e6
        lines!(ax, times ./ year, m; color = c, linestyle = kind == :VE ? :solid : :dash,
               label = kind == :VE ? sheet_name(k) : nothing)
    end
end
for (kind, r) in runs
    out = (injected .- [sum(co2_per_sheet(r, st)) for st in r.states]) ./ 1e6
    lines!(ax, times ./ year, out; color = :black,
           linestyle = kind == :VE ? :solid : :dash,
           label = kind == :VE ? "left the domain" : nothing)
end
axislegend(ax; position = :lt)
show_figure(fig, "layered_traps_sheets")

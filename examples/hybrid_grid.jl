# A hybrid grid (roadmap Stage 12): fine 3D cells around the injector, where
# the CO2 is not yet in vertical equilibrium, VE columns everywhere else. The
# same sloping, layered aquifer with a capillary fringe as
# layered_fringe.jl, run three ways: pure VE, hybrid, and 3D (every cell
# fine) as the reference.
#
#   julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/hybrid_grid.jl
#
# Self-contained (the grid is built here). One fine-scale fluid
# (FineScaleFluid) serves the VE columns (LayeredFringe, its upscaling) and
# the fine cells alike. Where a VE column meets a fine cell, the flux compares
# both sides' phase pressures at a common depth (a "virtual cell"; design note
# analysis/stage12-transition-flux-design.md).
#
# Fine cells must resolve the plume vertically: a fine cell taller than the
# CO2 layer smears it, and its low CO2 relperm then holds the plume back
# (design note §5). Here the cells are 0.25 m, the fringe about 1.7 m.
#
# Plotting backend: GLMakie (a window per figure) by default;
# JUTULCO2LAB_BACKEND=CairoMakie writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, JutulDarcy, Printf

year = 365.25 * 86400
g = Jutul.gravity_constant

# ------------------------ 1. The aquifer and its rock -------------------------

# As layered_fringe.jl, at twice its vertical resolution: a cross-section 3 km
# long, 10 m thick, 100 m wide, 100 columns of 40 cells (0.25 m), the caprock
# rising 1 % with gentle undulations; ten rock layers, permeability growing
# upwards from 30 mD to 1 D.
nx, nz, L, H, W = 100, 40, 3000.0, 10.0, 100.0
top(x) = 1100.0 - 0.01 * x + 1.5 * sin(2π * x / 750)
function aquifer_mesh()
    coord = Float64[]
    for j in 0:1, i in 0:nx
        x = i * L / nx
        append!(coord, [x, j * W, 0.0, x, j * W, 2000.0])
    end
    zc = [top((ii ÷ 2) * L / nx) + (k - 1 + s) * H / nz
          for k in 1:nz for s in (0, 1) for jj in 1:2 for ii in 1:2nx]
    grid = Dict("cartDims" => [nx, 1, nz], "COORD" => coord, "ZCORN" => zc,
                "ACTNUM" => fill(true, nx, 1, nz))
    return JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
end
mesh = aquifer_mesh()
nc = number_of_cells(mesh)
column(c) = cell_ijk(mesh, c)[1]
layer(c) = cell_ijk(mesh, c)[3]
md = 9.869233e-16
geo_layer(c) = cld(layer(c) * 10, nz)
perm = [30md * (1000 / 30)^((10 - geo_layer(c)) / 9) for c in 1:nc]
poro = [0.15 + 0.15 * (10 - geo_layer(c)) / 9 for c in 1:nc]

# ----------------------- 2. The fluid at the fine scale -----------------------

# A capillary fringe (brine saturation of capillary pressure, Pa), quadratic
# CO2 relperm, residual CO2 by co2lab's endpoint scaling.
rw, rg = 0.2, 0.2
fringe_curve = p -> rw + (1 - rw) * (5e3 / (max(p, 0.0) + 5e3))^2
fluid_fine = FineScaleFluid(;
    rock_types = [RockCurves(; invPc3D = fringe_curve, kr3D = s -> min(s, 1 - rw)^2)],
    residual = (rw, rg), default_pvt = nothing, co2_rho_pvt = nothing,
    co2_rho_ref = 700.0, wat_rho_pvt = (4.3e-10, 1e7), wat_rho_ref = 1050.0,
    co2_mu_pvt = nothing, co2_mu_ref = 6e-5, wat_mu_pvt = nothing, wat_mu_ref = 8e-4)

# ------------------------------- 3. Three grids -------------------------------

# Which cells are fine: none (pure VE), the pillars of columns 3–13 around the
# injector in column 8 (`ve_except_pillars` gives `ve_columns`' `is_ve` for
# whole pillars), or all (3D). Each fine cell becomes a block of its own.
well_column, fine_pillars = 8, 3:13
grids = ["VE" => fill(true, nc),
         "hybrid" => ve_except_pillars(mesh, (i, j) -> i in fine_pillars),
         "3D" => fill(false, nc)]

# --------------------------------- 4. One run ---------------------------------

# Injection for 2 years (about 20 % of the pore volume), then 8 years of
# migration; the up-dip end open (its column a brine reservoir, pore volume
# ×1e6). The injector: the VE column's, or the lowest 2 m of its fine cells.
rate = 0.2 * L * W * H * 0.2 / (2year)
dt = [fill(year / 12, 24); fill(year / 4, 32)]
injecting = [fill(true, 24); fill(false, 32)]

function run(is_ve)
    cols, _ = ve_columns(mesh; is_ve)
    ve = ve_domain(mesh, cols; perm_x = perm, poro)
    fluid = co2lab_ve_fluid(fluid_fine, mesh, cols; poro, perm_x = perm)
    blocks(i) = unique(cols.partition[c] for c in 1:nc if column(c) == i)
    wb = blocks(well_column)
    if length(wb) == 1
        well = setup_ve_well(ve, only(wb); name = :inj)
    else
        # Peaceman's index for each fine cell; a small wellbore volume (it
        # starts full of CO2)
        # (the lowest 2 m: nz/5 cells)
        perf = wb[end-nz÷5+1:end]
        cell(b) = cols.cells[cols.pos[b]]
        WI = [JutulDarcy.compute_peaceman_index(JutulDarcy.cell_dims(mesh, cell(b)),
                                                perm[cell(b)], 0.1) for b in perf]
        well = setup_well(ve.domain, perf; name = :inj, WI, simple_well = true,
                          volume_multiplier = 1e-3)
    end
    model = setup_co2lab_model(ve.domain; wells = [well], fluid)
    inject = co2_injector_control(TotalRateTarget(rate), fluid)
    on = setup_reservoir_forces(model; control = Dict(:inj => inject))
    off = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()))
    prm = setup_parameters(model)
    prm[:Reservoir][:StaticFluidVolume][blocks(nx)] .*= 1e6
    p0 = 1e7 .+ 1050.0 * g .* (ve.domain[:reference_depth] .- 1000.0)
    case = JutulCase(model, dt, [i ? on : off for i in injecting];
                     state0 = setup_co2lab_state(model; pressure = p0, forces = on),
                     parameters = prm)
    sim, cfg = setup_co2lab_simulator(case; info_level = -1)
    seconds = @elapsed res = simulate_reservoir(case; simulator = sim, config = cfg)
    its = sum(length(m[:steps]) for r in res.result.reports for m in r[:ministeps])
    return (ve = ve, columns = cols, fluid = fluid, states = res.states, its = its,
            seconds = seconds, blocks = number_of_columns(cols))
end

# (each run twice: the first compiles, the second is timed)
runs = [name => (run(is_ve); run(is_ve)) for (name, is_ve) in grids]

# --------------------------------- 5. Results ---------------------------------

# The CO2 saturation on the 3D cells (VE blocks reconstructed from their
# heights, fine cells their own), as a layer × column matrix
function saturation(r, k)
    s = ve_reconstruct(r.states[k], r.ve, r.columns, r.fluid).sG
    S = zeros(nz, nx)
    foreach(c -> S[layer(c), column(c)] = s[c], 1:nc)
    return S
end
# the L1 difference to 3D relative to the CO2 in 3D, by pore volume (all
# cells have the same volume); near the well (columns 1–15) and overall
pv = zeros(nz, nx)
foreach(c -> pv[layer(c), column(c)] = poro[c], 1:nc)
reference = last(runs[end])
function l1(r, k, cols = 1:nx)
    S, S3 = saturation(r, k)[:, cols], saturation(reference, k)[:, cols]
    return sum(pv[:, cols] .* abs.(S .- S3)) / sum(pv[:, cols] .* S3)
end
@printf("%-8s %7s %20s %20s %10s %8s\n", "", "blocks", "L1 near well (inj.)",
        "L1 all (inj. / end)", "Newton its", "time")
for (name, r) in runs
    @printf("%-8s %7d %20.3f %13.3f / %.3f %10d %7.2fs\n", name, r.blocks,
            l1(r, 24, 1:15), l1(r, 24), l1(r, length(dt)), r.its, r.seconds)
end

# --------------------------------- 6. Figures ---------------------------------

# Near the well at the end of injection, and the whole aquifer at the end:
# the hybrid's fine region (between the dashed lines) carries the 3D
# structure of the injected plume, its VE columns the VE one
dx = L / nx
fig = Figure(size = (1200, 750))
for (row, (name, r)) in enumerate(runs)
    for (col, (k, cols, when)) in enumerate(((24, 1:20, "2 years, near the well"),
                                             (length(dt), 1:nx, "10 years")))
        ax = Axis(fig[row, col]; title = "$name — $when", xlabel = "x (m)",
                  ylabel = "layer", yreversed = true)
        heatmap!(ax, ((cols) .- 0.5) .* dx, 1:nz, saturation(r, k)[:, cols]';
                 colorrange = (0, 0.8))
        edges = [first(fine_pillars) - 1, last(fine_pillars)] .* dx
        name == "hybrid" && vlines!(ax, edges; color = :white, linestyle = :dash)
    end
end
Colorbar(fig[:, 3]; colorrange = (0, 0.8), label = "CO2 saturation")
show_figure(fig, "hybrid_grid_plumes")

# The same panels on the grid itself: each cell as its own quadrilateral (its
# corners on the face y = 0), at its true depth. The layer-by-column panels
# above draw every cell as an equal rectangle, so the caprock's undulations
# show there as bumps in the saturation; here they are the grid's.
function cell_outline(c)
    bf = mesh.boundary_faces
    nodes = Set(n for f in bf.cells_to_faces[c] for n in bf.faces_to_nodes[f])
    pts = [Point2f(p[1], p[3]) for p in mesh.node_points[collect(nodes)] if p[2] == 0]
    # corners in order around the cell's centre
    m = sum(pts) / length(pts)
    return Makie.Polygon(sort(pts; by = p -> atan(p[2] - m[2], p[1] - m[1])))
end
outlines = cell_outline.(1:nc)
fig = Figure(size = (1200, 750))
for (row, (name, r)) in enumerate(runs)
    for (col, (k, xmax, when)) in enumerate(((24, 20dx, "2 years, near the well"),
                                             (length(dt), L, "10 years")))
        local ax = Axis(fig[row, col]; title = "$name — $when", xlabel = "x (m)",
                        ylabel = "depth (m)", yreversed = true)
        sG = ve_reconstruct(r.states[k], r.ve, r.columns, r.fluid).sG
        shown = findall(c -> column(c) * dx <= xmax, 1:nc)
        poly!(ax, outlines[shown]; color = sG[shown], colorrange = (0, 0.8))
        edges = [first(fine_pillars) - 1, last(fine_pillars)] .* dx
        name == "hybrid" && vlines!(ax, edges; color = :white, linestyle = :dash)
    end
end
Colorbar(fig[:, 3]; colorrange = (0, 0.8), label = "CO2 saturation")
show_figure(fig, "hybrid_grid_plumes_grid")

# The hybrid grid itself: its blocks outlined (black: the faces between cells
# of different blocks) on the cells' true geometry, the VE columns (blue) and
# the fine cells around the injector (orange), each fine cell a block of its
# own; near the well and over the whole aquifer
hybrid = last(runs[2])
# a face's two corners on the section y = 0 (none for faces lying in it)
face_edge(nodes) = [Point2f(p[1], p[3]) for p in mesh.node_points[collect(nodes)]
                    if p[2] == 0]
block_lines = Point2f[]
for f in 1:number_of_faces(mesh)
    l, r = mesh.faces.neighbors[f]
    hybrid.columns.partition[l] == hybrid.columns.partition[r] && continue
    e = face_edge(mesh.faces.faces_to_nodes[f])
    length(e) == 2 && append!(block_lines, e)
end
for f in eachindex(mesh.boundary_faces.neighbors)
    e = face_edge(mesh.boundary_faces.faces_to_nodes[f])
    length(e) == 2 && append!(block_lines, e)
end
is_fine = [!hybrid.columns.is_ve[hybrid.columns.partition[c]] for c in 1:nc]
fig = Figure(size = (1200, 650))
for (row, (xmax, title)) in enumerate(((20dx, "Hybrid grid near the well: VE columns " *
                                        "(blue) and fine cells (orange)"),
                                       (L, "The whole aquifer")))
    local ax = Axis(fig[row, 1]; title, xlabel = "x (m)", ylabel = "depth (m)",
                    yreversed = true)
    shown = findall(c -> column(c) * dx <= xmax, 1:nc)
    poly!(ax, outlines[shown]; color = [is_fine[c] ? :orange : :lightsteelblue
                                        for c in shown])
    keep = [k for k in 1:2:length(block_lines)
            if max(block_lines[k][1], block_lines[k+1][1]) <= xmax + 1e-6]
    linesegments!(ax, block_lines[[i for k in keep for i in (k, k + 1)]];
                  color = :black, linewidth = 0.5)
end
show_figure(fig, "hybrid_grid_blocks")

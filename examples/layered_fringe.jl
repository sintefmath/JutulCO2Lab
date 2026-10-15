# The capillary fringe in layered rock (LayeredFringe): CO2
# injected into a sloping aquifer whose permeability changes from layer to
# layer, and whose rock holds CO2 in a capillary fringe above the brine.
# co2lab can model the fringe or the layers, but not both at once; here the
# combined model runs next to the two single-effect ones and a 3D reference.
#
#   julia --project=examples examples/layered_fringe.jl
#
# Self-contained (the grid is built here). The fluid is described once at the
# fine scale (FineScaleFluid); the VE fluid is built from it, and so is the 3D
# reference (every cell fine, JutulCO2Lab's hybrid machinery).
#
# Plotting backend: GLMakie (a window per figure) by default;
# JUTULCO2LAB_BACKEND=CairoMakie writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, JutulDarcy, Printf

year = 365.25 * 86400
g = Jutul.gravity_constant

# ------------------------------- 1. The aquifer -------------------------------

# A 2D cross-section, 3 km long and 10 m thick, 100 columns of 20 cells (0.5
# m): the caprock rises 1 % towards x = 3 km with gentle undulations (small
# traps). Built as a corner-point grid, 100 m wide.
nx, nz, L, H, W = 100, 20, 3000.0, 10.0, 100.0
top(x) = 1100.0 - 0.01 * x + 1.5 * sin(2π * x / 750)
function aquifer_mesh()
    coord = Float64[]
    for j in 0:1, i in 0:nx
        x = i * L / nx
        append!(coord, [x, j * W, 0.0, x, j * W, 2000.0])
    end
    # each column's corners on its pillars: depth top(x) + (k − 1) H/nz
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

# ------------------------------ 2. Layered rock -------------------------------

# Ten geological layers of 1 m; permeability growing upwards (a prograding
# sequence: 30 mD at the bottom, 1 D at the top), porosity with it.
geo_layer(c) = cld(layer(c) * 10, nz)
md = 9.869233e-16 # millidarcy
perm = [30md * (1000 / 30)^((10 - geo_layer(c)) / 9) for c in 1:nc]
poro = [0.15 + 0.15 * (10 - geo_layer(c)) / 9 for c in 1:nc]

# ----------------------- 3. The fluid at the fine scale -----------------------

# Brine saturation of capillary pressure (Pa): from 1 at 0 down to the
# residual r_w, half-way at 5 kPa — a fringe of about 1.7 m for Δρ = 350
# kg/m³. CO2 relperm quadratic up to its endpoint, brine linear. Residual
# CO2 by co2lab's endpoint scaling. Constant densities and viscosities, the
# brine slightly compressible.
rw, rg = 0.2, 0.2
fringe_curve = p -> rw + (1 - rw) * (5e3 / (max(p, 0.0) + 5e3))^2
krg = s -> min(s, 1 - rw)^2
fluid_fine = FineScaleFluid(; rock_types = [RockCurves(; invPc3D = fringe_curve,
                                                       kr3D = krg)],
                            residual = (rw, rg), default_pvt = nothing,
                            co2_rho_pvt = nothing, co2_rho_ref = 700.0,
                            wat_rho_pvt = (4.3e-10, 1e7), wat_rho_ref = 1050.0,
                            co2_mu_pvt = nothing, co2_mu_ref = 6e-5,
                            wat_mu_pvt = nothing, wat_mu_ref = 8e-4)
# the same physics for co2lab's single-effect VE models
fluid_kw = (default_pvt = nothing, residual = (rw, rg), co2_rho_pvt = nothing,
            co2_rho_ref = 700.0, wat_rho_pvt = (4.3e-10, 1e7), wat_rho_ref = 1050.0,
            co2_mu_pvt = nothing, co2_mu_ref = 6e-5, wat_mu_pvt = nothing,
            wat_mu_ref = 8e-4)

# --------------------------------- 4. One run ---------------------------------

# Injection near the deep end (column 8) for 2 years, about 20 % of the pore
# volume, then 8 years of migration up-dip. The up-dip end is open: its
# column is a brine reservoir (pore volume ×1e6), so that CO2 leaves through
# an ordinary face, the same in VE and 3D.
rate = 0.2 * L * W * H * 0.2 / (2year)
dt = [fill(year / 12, 24); fill(year / 4, 32)]
injecting = [fill(true, 24); fill(false, 32)]
well_column = 8

# A run with the columns `cols` (VE blocks, or every cell fine) and VE fluid
# `fluid`: the domain, columns, fluid, states and Newton iterations
function run(cols, fluid)
    ve = ve_domain(mesh, cols; perm_x = perm, poro)
    blocks(i) = unique(cols.partition[c] for c in 1:nc if column(c) == i)
    # the injector: the VE column, or the lowest 2 m of fine cells with
    # Peaceman's index each (a small wellbore volume: it starts full of CO2)
    wb = blocks(well_column)
    if length(wb) == 1
        well = setup_ve_well(ve, only(wb); name = :inj)
    else
        perf = wb[end-3:end]
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
    res = simulate_reservoir(case; simulator = sim, config = cfg)
    its = sum(length(m[:steps]) for r in res.result.reports for m in r[:ministeps])
    return (ve = ve, columns = cols, fluid = fluid, states = res.states, its = its)
end

# ------------------------------- 5. Four models -------------------------------

columns, _ = ve_columns(mesh)
fine_columns, _ = ve_columns(mesh; is_ve = fill(false, nc))
models = [
    # the combined model, built from the fine-scale fluid
    "fringe + layers (LayeredFringe)" =>
        () -> run(columns, co2lab_ve_fluid(fluid_fine, mesh, columns; poro,
                                           perm_x = perm)),
    # co2lab's fringe: the column's rock averaged, one fringe for it
    "fringe only (CapillaryFringeTable)" =>
        () -> run(columns, co2lab_ve_fluid(; vertical = CapillaryFringeTable(;
                                                invPc3D = fringe_curve, kr3D = krg),
                                           fluid_kw...)),
    # co2lab's layered sharp interface: no capillarity
    "layers only (SharpInterfaceIntegrated)" =>
        () -> run(columns, co2lab_ve_fluid(; vertical = SharpInterfaceIntegrated(
                                               mesh, columns; poro, perm_x = perm),
                                           krmax = (1.0, krg(1 - rw)), fluid_kw...)),
    # the reference: every cell fine (3D), the same fine-scale fluid
    "3D (every cell fine)" =>
        () -> run(fine_columns, co2lab_ve_fluid(fluid_fine, mesh, fine_columns;
                                                poro, perm_x = perm))]
runs = [name => f() for (name, f) in models]

# --------------------------------- 6. Results ---------------------------------

# Each run's CO2 saturation on the 3D cells (VE: reconstructed from the
# heights), as a layer × column matrix
function saturation(r, k)
    s = ve_reconstruct(r.states[k], r.ve, r.columns, r.fluid).sG
    S = zeros(nz, nx)
    foreach(c -> S[layer(c), column(c)] = s[c], 1:nc)
    return S
end
# pore volume per cell (all cells have the same volume)
pv = zeros(nz, nx)
foreach(c -> pv[layer(c), column(c)] = poro[c], 1:nc)
reference = last(runs[end])
@printf("%-40s %16s %16s %10s\n", "", "L1 vs 3D (inj.)", "L1 vs 3D (end)",
        "Newton its")
# the L1 difference to 3D relative to the CO2 in 3D, by pore volume
function l1(r, k)
    S, S3 = saturation(r, k), saturation(reference, k)
    return sum(pv .* abs.(S .- S3)) / sum(pv .* S3)
end
for (name, r) in runs
    @printf("%-40s %16.3f %16.3f %10d\n", name, l1(r, 24), l1(r, length(dt)), r.its)
end

# --------------------------------- 7. Figures ---------------------------------

# The plume in the cross-section at the end of injection and at the end,
# for each model (depth below the caprock downwards, x along the aquifer)
fig = Figure(size = (1200, 900))
for (row, (name, r)) in enumerate(runs), (col, k) in enumerate((24, length(dt)))
    ax = Axis(fig[row, col]; title = "$name — $(k == 24 ? "2 years" : "10 years")",
              xlabel = "x (m)", ylabel = "layer", yreversed = true)
    heatmap!(ax, range(0, L; length = nx), 1:nz, saturation(r, k)';
             colorrange = (0, 0.8))
end
Colorbar(fig[:, 3]; colorrange = (0, 0.8), label = "CO2 saturation")
show_figure(fig, "layered_fringe_plumes")

# One column's saturation profile at the end of injection: the fringe's
# smooth tail, and the layers' jumps (high permeability holds more CO2 near
# the caprock)
fig = Figure(size = (700, 600))
ax = Axis(fig[1, 1]; title = "Column 20 at the end of injection",
          xlabel = "CO2 saturation", ylabel = "layer", yreversed = true)
for (name, r) in runs
    lines!(ax, saturation(r, 24)[:, 20], 1:nz; label = name)
end
axislegend(ax; position = :rb)
show_figure(fig, "layered_fringe_profile")

# The same plumes drawn on the grid itself: each cell as its own quadrilateral
# (its corners on the face y = 0), so the undulating caprock shows as it is
function cell_outline(c)
    bf = mesh.boundary_faces
    nodes = Set(n for f in bf.cells_to_faces[c] for n in bf.faces_to_nodes[f])
    pts = [Point2(p[1], p[3]) for p in mesh.node_points[collect(nodes)] if p[2] == 0]
    # corners in order around the cell's centre
    m = sum(pts) / length(pts)
    return Makie.Polygon(sort(pts; by = p -> atan(p[2] - m[2], p[1] - m[1])))
end
outlines = cell_outline.(1:nc)
fig = Figure(size = (1200, 900))
for (row, (name, r)) in enumerate(runs), (col, k) in enumerate((24, length(dt)))
    local ax = Axis(fig[row, col]; xlabel = "x (m)", ylabel = "depth (m)",
                    title = "$name — $(k == 24 ? "2 years" : "10 years")",
                    yreversed = true)
    sG = ve_reconstruct(r.states[k], r.ve, r.columns, r.fluid).sG
    poly!(ax, outlines; color = sG, colorrange = (0, 0.8))
end
Colorbar(fig[:, 3]; colorrange = (0, 0.8), label = "CO2 saturation")
show_figure(fig, "layered_fringe_plumes_grid")

# Every report step of every model in Jutul's interactive viewer (GLMakie
# only): the menu at the bottom right picks the model, the slider the step
if backend == :GLMakie
    sat = [name => [ve_reconstruct(s, r.ve, r.columns, r.fluid).sG for s in r.states]
           for (name, r) in runs]
    states = [Dict(Symbol(name) => s[k] for (name, s) in sat) for k in eachindex(dt)]
    fig = plot_interactive(mesh, states; new_window = false, axis_view = :xz,
                           aspect = (6, 0.1, 1))
    show_figure(fig, "layered_fringe_viewer")
end

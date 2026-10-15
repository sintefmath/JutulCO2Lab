# Grid models from the CO2 Storage Atlas (co2lab's modelsFromAtlas.m): a
# formation as a volumetric grid, all North Sea formations in 3D, and their
# volumes.
#
#   julia --project=examples examples/atlas/models_from_atlas.jl
#
# Plotting backend: GLMakie (a window per figure; close it to continue) by
# default; JUTULCO2LAB_BACKEND=CairoMakie writes the figures to files instead.

backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using JutulCO2Lab, Jutul, Printf

# ------------------------------ 1. Utsira in 3D -------------------------------

# The top-surface and thickness maps are interpolated onto one lattice, and
# every cell where both are defined (and at least 1 mm thick) becomes a cell
# of a corner-point grid; `layers` would split each into several. Utsira at
# full resolution has about 100 000 cells; taking every second map value in
# each direction leaves about 25 000.
m = atlas_mesh("Utsirafm"; coarsening = 2, layers = 1)
println("Utsira, coarsening 2: ", number_of_cells(m.mesh), " cells")

# Coarsening can split a formation into pieces; atlas_mesh keeps the largest.
# The grid, coloured by cell volume: the folds of the top surface are the
# targets for structural trapping of migrating CO2.
volumes = tpfv_geometry(m.mesh).volumes
fig = Figure(size = (1100, 900))
ax = Axis3(fig[1, 1]; title = "Utsira: cell volumes (m³)", zreversed = true,
           aspect = (1, 3.5, 0.6), azimuth = -0.4, elevation = 0.9)
p = plot_cell_data!(ax, m.mesh, volumes)
Colorbar(fig[1, 2], p)
show_figure(fig, "utsira_cell_volumes")

# ----------------------- 2. All North Sea formations --------------------------

# Only formations with both a top and a thickness map (their own, or the
# Jurassic top) can be built as volumetric grids.
names = atlas_grid_formations(:north_sea)

# Each formation's grid, at every second map value, in a panel of its own,
# coloured by depth.
fig = Figure(size = (1600, 1600))
for (k, name) in enumerate(names)
    g = atlas_mesh(name; coarsening = 2).mesh
    depth = tpfv_geometry(g).cell_centroids[3, :]
    panel = Axis3(fig[fldmod1(k, 4)...]; title = name, zreversed = true,
                  aspect = (1, 1, 0.3), azimuth = -1.0, elevation = 0.7)
    hidedecorations!(panel)
    plot_cell_data!(panel, g, depth; colormap = :turbo)
end
show_figure(fig, "north_sea_grids")

# ----------------------------- 3. Their volumes -------------------------------

# Full resolution: number of cells, depth range of the caprock and bulk
# volume. co2lab adds the volume in structural traps, from the trap analysis
# (see the examples in trapping/).
println("\nNorth Sea formations, full resolution:")
println("  name          |  cells | top min | top max | bulk volume (m³)")
println("  --------------|--------|---------|---------|-----------------")
for name in names
    g = atlas_mesh(name)
    z_top = ve_column_geometry(g.mesh, g.columns).z_top
    volume = sum(tpfv_geometry(g.mesh).volumes)
    @printf("  %-13s | %6d | %7.0f | %7.0f | %10.3e\n", name,
            number_of_cells(g.mesh), extrema(z_top)..., volume)
end

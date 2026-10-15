# The CO2 Storage Atlas of the Norwegian North Sea (co2lab's showCO2atlas.m):
# the raw maps, the formations' outlines on a map of Norway with the shelf's
# wells, and the formations in 3D.
#
#   julia --project=examples examples/atlas/show_co2_atlas.jl
#
# Plotting backend: GLMakie (interactive windows) by default; set
# JUTULCO2LAB_BACKEND=CairoMakie to write the figures to files instead (no
# display needed). Data: the Atlas data sets in atlas_data_directory().

backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using JutulCO2Lab, Jutul, Printf

# ------------------------------- 1. The raw data ------------------------------

# Each formation comes as maps on regular lattices: depth of the top surface
# and/or thickness, at 200-1000 m resolution, as NPD published them.
maps = read_atlas_maps(atlas_formation_names(:north_sea))
println("North Sea maps:")
for m in maps
    @printf("  %-14s %-9s %4.0f m, %4d x %-4d values\n", m.name, m.variant,
            m.meta.cellsize, m.meta.dims...)
end

# Some formations have only a thickness map; their top is taken from the
# Jurassic maps. Those with both kinds of map can be built as grids:
names = atlas_grid_formations(:north_sea)
println("\nformations with a grid: ", join(names, ", "))

# ------------------------------ 2. Utsira, raw ------------------------------

# The maps are usable on their own. Note that the top and thickness maps of a
# formation need not cover the same area or share a lattice. Each is drawn as
# a surface with its values as height (depth downwards), seen from above as
# co2lab's surf(...) with view(0, 90) — rotatable with the mouse in GLMakie.
fig = Figure(size = (900, 900))
for (k, m) in enumerate(filter(m -> m.name == "Utsirafm", maps))
    ax = Axis3(fig[1, 2k - 1]; title = "Utsira $(m.variant) (m)",
               zreversed = m.variant == :top, aspect = (1, 3, 0.3),
               azimuth = -pi / 2, elevation = pi / 2, xticks = LinearTicks(3),
               zticklabelsvisible = false, zlabelvisible = false)
    plot_atlas_map!(ax, m)
    Colorbar(fig[1, 2k], ax.scene.plots[end])
end
resize_to_layout!(fig)
show_figure(fig, "utsira_raw_maps")

# ----------------------- 3. The formations on the map -------------------------

# Each formation as a grid, 5 x 5 map values per cell, drawn by its outline at
# its depth, with Norway's coastline and the wells of the shelf (at sea level)
# for orientation; seen from above, as co2lab's view(2).
fig = Figure(size = (1100, 1000))
ax = Axis3(fig[1, 1]; title = "North Sea formations", zreversed = true,
           aspect = (1, 1.25, 0.2), azimuth = -pi / 2, elevation = pi / 2,
           limits = ((2.5e5, 8.5e5), (6.2e6, 6.95e6), nothing),
           xticks = LinearTicks(3), zticklabelsvisible = false, zlabelvisible = false)
colors = Makie.to_colormap(:tab20)                  # 20 distinct colours
for (k, name) in enumerate(names)
    m = atlas_mesh(name; coarsening = 5)
    surface = ve_top_surface(m.mesh, m.columns)
    plot_ve_outline!(ax, surface; color = colors[k], linewidth = 2, label = name)
end
plot_atlas_context!(ax)
Legend(fig[1, 2], ax)
show_figure(fig, "north_sea_formations")

# ------------------------------ 4. The same in 3D -----------------------------

# The formations at their depths (3 x 3 map values per cell), each in one
# colour, with Jutul's own mesh plotting.
fig = Figure(size = (1300, 900))
ax = Axis3(fig[1, 1]; title = "North Sea formations in 3D", zreversed = true,
           aspect = (1, 1.4, 0.4), azimuth = 1.2, elevation = 0.5)
for (k, name) in enumerate(names)
    m = atlas_mesh(name; coarsening = 3)
    plot_cell_data!(ax, m.mesh, fill(k, number_of_cells(m.mesh));
                    colormap = :hsv, colorrange = (1, length(names)))
end
show_figure(fig, "north_sea_formations_3d")

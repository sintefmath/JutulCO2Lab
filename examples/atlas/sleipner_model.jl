# The Sleipner benchmark model (co2lab's makeSleipnerVEmodel): layer 9 of the
# Utsira formation around the Sleipner injection site, from the IEAGHG
# benchmark, as a VE model in global coordinates — shown in 3D and on the map
# of the Atlas's Utsira formation.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/atlas/sleipner_model.jl
#
# The benchmark data are no longer downloadable: a local copy must be in
# sleipner_data_directory(). Plotting backend: GLMakie (a window per figure;
# close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie writes the
# figures to files instead.

using JutulCO2Lab

# ------------------------------- 0. The data ----------------------------------

# Stop with a message, rather than fail later, if the data are not here.
if !isfile(joinpath(sleipner_data_directory(), "M9X1.grdecl"))
    println(stderr, "The Sleipner benchmark data are not in ", sleipner_data_directory(),
            ".\nThe IEAGHG benchmark model is no longer publicly downloadable; " *
            "place a copy of its files (M9X1.grdecl, M9X1_*.inc) there.")
    exit(1)
end

backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, Printf

# ------------------------------ 1. The model ----------------------------------

# The benchmark grid (65 x 119 x 43), its shale cells (below 200 mD) removed,
# placed in UTM zone 31 by the file's MAPAXES; columns and the VE domain.
s = sleipner_ve_domain()
geo, rock = s.ve.geometry, s.ve.rock
@printf("%d cells in %d columns; caprock %.0f-%.0f m, height %.1f-%.1f m\n",
        number_of_cells(s.mesh), number_of_columns(s.columns), extrema(geo.z_top)...,
        extrema(geo.H)...)
@printf("column averages: permeability %.2f-%.2f D, porosity %.3f-%.3f\n",
        (extrema(rock.perm) ./ 9.869233e-13)..., extrema(rock.poro)...)
println("open boundary: ", length(s.sides), " lateral column sides")

# ------------------------------- 2. In 3D -------------------------------------

# The grid coloured by depth (as makeSleipnerVEmodel shows it).
depth = tpfv_geometry(s.mesh).cell_centroids[3, :]
fig = Figure(size = (1100, 900))
ax = Axis3(fig[1, 1]; title = "Sleipner benchmark, layer 9: depth (m)", zreversed = true,
           aspect = (1, 1.8, 0.5), azimuth = -0.9, elevation = 0.6)
p = plot_cell_data!(ax, s.mesh, depth; colormap = :turbo)
Colorbar(fig[1, 2], p)
show_figure(fig, "sleipner_3d")

# ------------------------ 3. On the Atlas's Utsira map ------------------------

# In global coordinates the benchmark can be set against the Atlas: the
# Atlas's Utsira top depth at the benchmark's column centres, against the
# benchmark's caprock.
surface = ve_top_surface(s.mesh, s.columns)
top = only(filter(m -> m.variant == :top, read_atlas_maps(["Utsirafm"])))
X, Y = JutulCO2Lab._map_positions(top.meta)
atlas_depth = [JutulCO2Lab._interp2(X, Y, permutedims(top.data), x, y)
               for (x, y) in zip(surface.xc, surface.yc)]
d = geo.z_top .- atlas_depth
@printf("benchmark caprock minus Atlas top: mean %.1f m, range %.1f to %.1f m\n",
        sum(d) / length(d), extrema(d)...)

# The Atlas's Utsira around Sleipner, coloured by top depth, with the
# benchmark's outline; next to it the benchmark's own caprock.
utsira = atlas_mesh("Utsirafm")
fig = Figure(size = (1400, 700))
ax1 = Axis(fig[1, 1]; title = "Atlas: Utsira top depth (m), benchmark outlined",
           aspect = DataAspect(), limits = ((4.25e5, 4.55e5), (6.46e6, 6.485e6)))
us = ve_top_surface(utsira.mesh, utsira.columns)
plot_ve_map!(ax1, us; background = ve_column_geometry(utsira.mesh, utsira.columns).z_top,
             depth_contours = 0, outline = false)
plot_ve_outline!(ax1, surface; color = :red)
ax2 = Axis(fig[1, 2]; title = "benchmark caprock depth (m)", aspect = DataAspect())
plot_ve_map!(ax2, surface; background = geo.z_top, depth_contours = 10)
for (k, ax) in enumerate((ax1, ax2))
    Colorbar(fig[2, k], ax.scene.plots[findfirst(p -> p isa Poly, ax.scene.plots)];
             vertical = false)
end
show_figure(fig, "sleipner_on_utsira")

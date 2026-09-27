# co2lab-spillpoint's subscaleResolution: how much trapping volume the CO2
# Storage Atlas misses because its grids cannot resolve small-scale caprock
# relief. The Sleipner benchmark model (about 50 m cells) resolves traps the
# Atlas's Utsira grid (500 m) cannot; its trap volume per area, scaled to
# Utsira's area, estimates the volume unresolved in all of Utsira — first as
# it is, then with the large-scale trend of the caprock removed.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/trapping/subscale_resolution.jl
#
# The Sleipner benchmark data are no longer downloadable: a local copy must be
# in sleipner_data_directory(). Plotting backend: GLMakie (a window per
# figure, 3D views rotatable with the mouse; close it to continue) by default;
# JUTULCO2LAB_BACKEND=CairoMakie writes the figures to files instead.

using JutulCO2Lab

# Stop with a message, rather than fail later, if the data are not here.
if !isfile(joinpath(sleipner_data_directory(), "M9X1.grdecl"))
    println(stderr, "The Sleipner benchmark data are not in ", sleipner_data_directory(),
            ".\nThe IEAGHG benchmark model is no longer publicly downloadable; " *
            "place a copy of its files (M9X1.grdecl, M9X1_*.inc) there.")
    exit(1)
end

backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using JutulDarcy, Printf, Statistics

# ----------------------------- 1. The two models ------------------------------

# The whole Sleipner benchmark grid (shale cells included, as co2lab reads it
# here) in global coordinates, and the Atlas's Utsira formation.
grid = sleipner_grid_section(; min_perm = 0.0)
mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
sleipner = trap_surface(mesh, first(ve_columns(mesh)))
utsira_grid = atlas_mesh("Utsirafm")
utsira = trap_surface(utsira_grid.mesh, utsira_grid.columns)

# Cell-based traps on both.
res_sleipner = trap_analysis(sleipner; method = :cell)
res_utsira = trap_analysis(utsira; method = :cell)
@printf("Sleipner: %d columns, %d traps; Utsira: %d columns, %d traps\n",
        length(sleipner.zc), maximum(res_sleipner.traps), length(utsira.zc),
        maximum(res_utsira.traps))

# ------------------------- 2. Unresolved trap volume --------------------------

# Sleipner's trap volume per area of caprock, applied to all of Utsira.
area_utsira = sum(utsira.area)
trapvol_utsira = sum(trap_volumes(utsira, res_utsira))
function unresolved(s, res)
    per_area = sum(trap_volumes(s, res)) / sum(sleipner.area)
    lost = area_utsira * per_area
    @printf("  %.4g litres of trap volume per m2 not resolved by the Atlas grid\n",
            1000 * per_area)
    @printf("  in all of Utsira %.4g litres, %.2f %% of its resolved trap volume\n",
            1000 * lost, 100 * lost / trapvol_utsira)
    return lost
end
println("Sleipner as it is:")
unresolved(sleipner, res_sleipner)

# The Utsira columns under Sleipner's extent (red) and Sleipner's traps.
region = (minimum(sleipner.x) .< utsira.cx .< maximum(sleipner.x)) .&
         (minimum(sleipner.y) .< utsira.cy .< maximum(sleipner.y))
fig = Figure(size = (1100, 900))
ax = ve_axis3(fig[1, 1], sleipner; view = (-40, 50), exaggeration = 5,
              title = "Utsira under Sleipner (red) and Sleipner's traps")
plot_ve_surface!(ax, (; utsira..., polygons = utsira.polygons[region]);
                 color = (:red, 0.4), shading = true, strokewidth = 0)
plot_traps!(ax, sleipner, res_sleipner; rivers = false)
show_figure(fig, "sleipner_on_utsira")

# ---------------------- 3. The large-scale trend removed ----------------------

# Part of Sleipner's relief is Utsira's own large-scale trend, which the
# Atlas resolves. It is removed by subtracting the Atlas's Utsira top depth
# (interpolated bilinearly; co2lab interpolates the Utsira columns' depths on
# a Delaunay triangulation) and adding back its mean over Sleipner's extent.
top = only(filter(m -> m.variant == :top, read_atlas_maps(["Utsirafm"])))
level = mean(utsira.zc[region])
trend(x, y) = atlas_map_values(top, x, y) .- level
detrended = (; sleipner..., zc = sleipner.zc .- trend(sleipner.cx, sleipner.cy),
             z = sleipner.z .- trend(sleipner.x, sleipner.y))
res_detrended = trap_analysis(detrended; method = :cell)
println("Sleipner with the trend removed ($(maximum(res_detrended.traps)) traps):")
lost = unresolved(detrended, res_detrended)

# Sleipner's traps before and after, next to the share of trap volume in
# local and in large-scale relief.
fig = Figure(size = (1300, 750))
for (k, (name, s, r)) in enumerate((("as it is", sleipner, res_sleipner),
                                    ("trend removed", detrended, res_detrended)))
    a = ve_axis3(fig[1, k], s; view = (-50, 40), exaggeration = 5,
                 title = "Sleipner traps (red), $name")
    plot_ve_map!(a, s; background = s.zc, regions = r.traps .> 0, region_color = :red)
end
ax = Axis(fig[1, 3][1, 1]; title = "Utsira trap volume", aspect = 1)
pie!(ax, [lost, trapvol_utsira]; color = [:orange, :steelblue])
hidedecorations!(ax); hidespines!(ax)
Legend(fig[1, 3][2, 1], [PolyElement(color = c) for c in (:orange, :steelblue)],
       ["local variations (unresolved)", "large-scale variations"])
resize_to_layout!(fig)
show_figure(fig, "sleipner_subscale")

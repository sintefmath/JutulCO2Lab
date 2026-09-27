# co2lab-spillpoint's coarseningEffects: how coarsening the geometry changes
# structural trapping, on the Johansen formation of the CO2 Storage Atlas at
# six resolutions — the full data set (500 m cells, about 80 000 columns) and
# coarsened by 2 to 6 (3000 m cells, about 2000 columns). Coarsening smooths
# the caprock; the small ridges and domes that divert CO2 and trap it are the
# first to go.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/trapping/coarsening_effects.jl
#
# Plotting backend: GLMakie (a window per figure, 3D views rotatable with the
# mouse; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Printf, Statistics

# ----------------------------- 1. Six resolutions -----------------------------

# Each resolution's caprock surface, its cell-based traps and their volumes.
N = 6
cases = map(1:N) do c
    m = atlas_mesh("Johansenfm"; coarsening = c)
    s = trap_surface(m.mesh, m.columns)
    res = trap_analysis(s; method = :cell)
    volumes = trap_volumes(s, res)
    @printf("coarsening %d (%4d m): %6d columns, depth std %.4f m, ", c, 500c,
            length(s.zc), std(s.zc))
    @printf("%4d traps, %.3e m3\n", length(volumes), sum(volumes))
    (surface = s, res = res, volumes = volumes)
end

# ------------------------------ 2. The subdomain ------------------------------

# A box of 5 km x 50 km, shown on the full-resolution caprock.
box = ((5.25e5, 5.30e5), (6.70e6, 6.75e6))
inbox(s) = (box[1][1] .< s.cx .< box[1][2]) .& (box[2][1] .< s.cy .< box[2][2])
s1 = cases[1].surface
fig = Figure(size = (900, 1000))
ax = ve_axis3(fig[1, 1], s1; view = (0, 90), exaggeration = 5,
              title = "Johansen: caprock depth, the subdomain shaded")
plot_ve_map!(ax, s1; background = s1.zc, regions = inbox(s1),
             region_color = (:black, 0.5))
show_figure(fig, "johansen_subdomain")

# ---------------- 3. Traps and trap volumes against resolution ----------------

# The subdomain at each resolution, 1000 m apart in depth, with the traps
# filled to their spill depth (co2lab's flattenTraps) in solid colour: many
# small traps disappear as the resolution coarsens.
colours = [:blue, :red, :green, :cyan, :magenta, :yellow]
fig = Figure(size = (1500, 900))
# (co2lab's view(86, 12): azimuth 86° − 90° from x)
ax = Axis3(fig[1:2, 1]; title = "traps in the subdomain, coarsening 1 (top) to 6",
           zreversed = true, aspect = (1, 4, 3), azimuth = deg2rad(-4),
           elevation = deg2rad(12))
for (c, case) in enumerate(cases)
    s, sub = case.surface, inbox(case.surface)
    part = (; s..., polygons = s.polygons[sub])
    plot_ve_surface!(ax, part; color = (colours[c], 0.2), offset = 1000.0 * c,
                     strokecolor = (colours[c], 0.6), strokewidth = 0.5)
    flat = flatten_traps(s, case.res)
    traps = (; flat..., polygons = s.polygons[sub .& (case.res.traps .> 0)])
    plot_ve_surface!(ax, traps; color = colours[c], offset = 1000.0 * c,
                     strokewidth = 0)
end

# Next to it, the trap volumes of the whole formation. The total volume first
# grows (the largest traps grow a little at lower resolution), then shrinks as
# smaller traps vanish; the mean volume grows quickly as the small traps are
# smoothed away.
labels = (1:N, string.(500 .* (1:N)))
barplot(fig[1, 2], 1:N, [sum(c.volumes) for c in cases]; color = colours,
        axis = (title = "total trap volume (m³)", xticks = labels,
                xlabel = "cell size (m)"))
barplot(fig[2, 2], 1:N, [mean(c.volumes) for c in cases]; color = colours,
        axis = (title = "mean trap volume (m³)", xticks = labels,
                xlabel = "cell size (m)"))
colsize!(fig.layout, 1, Relative(0.6))
show_figure(fig, "johansen_coarsening")

# co2lab-spillpoint's optimizationExample: where to inject into the Utsira
# formation to fill as much structural trap volume as possible. Traps spill
# into each other along trees; injecting below a tree's root trap fills the
# whole tree as CO2 spills onwards. The trees are ranked by volume, and the
# best single injection point is found.
#
#   julia --project=examples examples/trapping/optimization_example.jl
#
# Plotting backend: GLMakie (a window per figure, 3D views rotatable with the
# mouse; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Printf

# ------------------- 1. The Utsira formation and its traps --------------------

# The Atlas's Utsira grid coarsened by 2 (1 km cells), cell-based traps.
m = atlas_mesh("Utsirafm"; coarsening = 2)
s = trap_surface(m.mesh, m.columns)
res = trap_analysis(s; method = :cell)
@printf("Utsira: %d columns, %d traps (average porosity %.2f)\n", length(s.zc),
        maximum(res.traps), m.petro.poro)

# ------------------------------- 2. Spill trees -------------------------------

# A tree: a root trap nothing spills into and every trap downstream of it.
# Trees overlap where they drain into the same traps; without overlap each
# trap counts in the largest tree that holds it (volumes in bulk m3).
trees = maximize_trapping(s, res)
trees_overlap = maximize_trapping(s, res; remove_overlap = false)
tv = [t.value for t in trees]
with_overlap = Dict(t.root => t.value for t in trees_overlap)
extra = [with_overlap[t.root] - t.value for t in trees]
@printf("%d trees; the largest %.3e m3, %.1f %% of all trap volume\n", length(trees),
        tv[1], 100 * tv[1] / sum(tv))

fig = Figure(size = (1100, 500))
ax = Axis(fig[1, 1]; yscale = log10, xlabel = "tree", ylabel = "tree volume (m³)")
barplot!(ax, 1:length(tv), tv; label = "without overlap")
barplot!(ax, 1:length(tv), tv .+ extra; fillto = tv, label = "overlap")
ylims!(ax, minimum(tv), 1.1 * maximum(tv .+ extra))
axislegend(ax; position = :rt)
show_figure(fig, "utsira_trees")

# ---------------------------- 3. Injection points -----------------------------

# For the ten largest trees, the deepest column of the root's spill region:
# CO2 injected there rises through the region into the root trap. The best
# single point (co2lab's findOptimalInjectionPoint) lies on the border of
# two spill regions, so that it fills both of their trees.
N = 10
roots = [t.root for t in trees[1:N]]
deepest(r) = (cols = findall(==(r), res.trap_regions); cols[argmax(s.zc[cols])])
points = deepest.(roots)
best = optimal_injection_point(s, res)
@printf("best single injection point: column %d, %.3e m3 (%.1f %% of the volume)\n",
        best.column, best.volume, 100 * best.volume / sum(tv))

# Each trap coloured by log10 of its tree's volume, the injection points with
# the share of the total trap volume they can reach.
tree_value = zeros(maximum(res.traps))
for t in reverse(trees)
    tree_value[t.traps] .= log10(t.value)
end
fig = Figure(size = (1100, 1000))
ax = ve_axis3(fig[1, 1], s; view = (64, 72), exaggeration = 40,
              title = "Utsira: spill trees and injection points")
plot_traps!(ax, s, res; trap_values = tree_value, rivers = false)
Colorbar(fig[1, 2]; limits = extrema(tree_value), label = "log₁₀ tree volume (m³)")
# the points a little above the caprock
xy(c) = Point3f(s.cx[c], s.cy[c], s.zc[c] - 20)
scatter!(ax, [xy(c) for c in points]; color = :red, strokecolor = :blue,
         strokewidth = 2, markersize = 12)
text!(ax, [xy(c) for c in points]; text = [@sprintf("%.1f%%", 100 * v / sum(tv))
                                          for v in tv[1:N]], offset = (6, 0))
scatter!(ax, [xy(best.column)]; color = :yellow, strokecolor = :black,
         strokewidth = 2, markersize = 16, marker = :star5)
best_share = 100 * best.volume / sum(tv)
text!(ax, [xy(best.column)]; text = @sprintf("best: %.1f%%", best_share),
      offset = (8, -12))
show_figure(fig, "utsira_injection_points")

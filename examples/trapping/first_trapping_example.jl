# co2lab-spillpoint's firstTrappingExample: the basic routines of the trapping
# analysis on a synthetic shoe-box — structural traps, the spill paths
# ("rivers") between them, spill regions, the trapping capacity, and the leaf
# traps, which nothing spills into.
#
#   julia --project=examples examples/trapping/first_trapping_example.jl
#
# Plotting backend: GLMakie (a window per figure, 3D views rotatable with the
# mouse; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead. It ends, as co2lab's, with the
# interactive trap viewer.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, Printf

# ---------------------- 1. The grid and its top surface -----------------------

# A sandbox of 10 km x 5 km x 50 m at 1000 m depth, dipping in x, its top
# perturbed by sin/cos terms that make domes; 100 x 100 cells, one layer.
Lx, Ly, H = 10000.0, 5000.0, 50.0
x, y = range(0, Lx; length = 101), range(0, Ly; length = 101)
bump(a, b) = a - 0.2 * sin(5pi * a) * sin(5pi * b^1.5) - 0.075 * sin(1.25pi * b) +
             0.15 * sin(a + b)
top = [1000 + H * bump(a / Lx, b / Ly) for a in x, b in y]
bottom = [1000 + H * (1 + a / Lx) for a in x, b in y]
grid = surface_mesh(x, y, top; bottom)

# The top surface as the trapping analysis sees it: column polygons, their
# caprock depths, and the nodes and edges of the surface.
s = trap_surface(grid.mesh, grid.columns)

# -------------------------- 2. Traps and spill paths --------------------------

# The cell-based method (co2lab's trapAnalysis(Gt, true)): traps from the
# columns' depths and adjacency. Each trap spills over its spill point, and CO2
# then follows a river to the next trap or out of the domain.
res = trap_analysis(s; method = :cell)
ntrap = maximum(res.traps)
@printf("%d traps; %d of %d columns in traps\n", ntrap, count(>(0), res.traps),
        length(res.traps))

# The 3D grid, and above it (lifted 15 m) the top surface coloured by caprock
# depth with the traps in white — rotatable with the mouse in GLMakie.
fig = Figure(size = (1100, 700))
ax = ve_axis3(fig[1, 1], s; view = (30, 25), exaggeration = 25,
              title = "grid, caprock depth and traps (white)")
plot_mesh!(ax, grid.mesh; color = :lightyellow, transparency = true, alpha = 0.3)
plot_mesh_edges!(ax, grid.mesh; color = (:black, 0.1))
plot_ve_map!(ax, (; s..., z = s.z .- 15); background = s.zc, colormap = :jet,
             regions = res.traps .> 0, region_color = :white, outline = false)
show_figure(fig, "shoebox_3d")

# -------------------------- 3. The trapping capacity --------------------------

# The pore volume of each trap (porosity 0.25), largest first, against the
# pore volume of the whole sandbox.
poro = fill(0.25, length(s.zc))
volumes = trap_volumes(s, res; poro)
total = sum(volumes)
pv = sum(s.area .* s.H .* poro)
@printf("total trapping capacity %.3e m3, %.2f %% of the pore volume %.2e m3\n\n",
        total, 100 * total / pv, pv)
order = sortperm(volumes; rev = true)
ncols = [count(==(t), res.traps) for t in order]
println("trap   | volume (m3) | columns")
for (t, n) in zip(order, ncols)
    @printf("%6d | %11.3e | %7d\n", t, volumes[t], n)
end
@printf("\nthe five largest traps hold %.2e m3, %.1f %% of the capacity\n",
        sum(volumes[order[1:5]]), 100 * sum(volumes[order[1:5]]) / total)

# Traps, rivers and spill regions, seen from above, with trap volumes and
# sizes.
fig = Figure(size = (1500, 650))
ax = ve_axis3(fig[1:2, 1], s; view = (90, 90), exaggeration = 25,
              title = "traps, spill paths and spill regions")
plot_traps!(ax, s, res; labels = true)
labels = (1:ntrap, string.(order))
barplot(fig[1, 2], 1:ntrap, volumes[order];
        axis = (title = "trap volume (m³)", xticks = labels))
barplot(fig[2, 2], 1:ntrap, ncols; axis = (title = "columns in trap", xticks = labels))
colsize!(fig.layout, 1, Relative(0.55))
show_figure(fig, "shoebox_capacity")

# ------------------------------- 4. Leaf traps --------------------------------

# Leaf traps: traps that nothing spills into. CO2 injected in a leaf's spill
# region fills the leaf first and then follows its rivers upward, until it
# reaches a trap that collects CO2 from several leaves or leaves the domain.
into = vec(sum(res.trap_adj; dims = 1))
leaves = findall(==(0), into)
leaf_volumes = trap_volumes(s, res; traps = leaves)
for (t, v) in zip(leaves, leaf_volumes)
    next = findall(res.trap_adj[t, :])
    @printf("leaf trap %2d: bulk volume %.3e m3, spills %s\n", t, v,
            isempty(next) ? "out of the domain" : "into trap(s) " * join(next, ", "))
end

# Each leaf's spill region and rivers in a colour of its own, the leaves white,
# the other traps dark grey.
fig = Figure(size = (1100, 650))
ax = ve_axis3(fig[1, 1], s; view = (120, 25), exaggeration = 25,
              title = "leaf traps (white), their spill regions and rivers")
leaf_regions = [t in leaves ? t : 0 for t in res.trap_regions]
drawn = (; res..., trap_regions = leaf_regions, traps = zeros(Int, length(res.traps)),
         rivers = [t in leaves ? res.rivers[t] : Vector{Int}[] for t in 1:ntrap])
plot_traps!(ax, s, drawn)
is_leaf = res.traps .∈ Ref(leaves)
plot_ve_map!(ax, s; regions = (res.traps .> 0) .& .!is_leaf,
             region_color = (:gray30, 1.0), depth_contours = 0, outline = false)
plot_ve_map!(ax, s; regions = is_leaf, region_color = :white,
             depth_contours = 0, outline = false)
show_figure(fig, "shoebox_leaves")

# ----------------------- 5. The interactive trap viewer -----------------------

# The shoe-box in the interactive viewer (interactive_traps, co2lab's
# interactiveTrapping), as co2lab opens it: cell-based traps, lit, spill
# regions shown, the path traps in grey, column 1898 selected. Click a column
# to see where CO2 injected there goes; see show_traps_interactively.jl.
v = interactive_traps(s; method = :cell, lighting = true, spill_regions = true,
                      path_colours = false, injection_column = 1898, view = (80, 20),
                      title = "Shoe-box: interactive trap viewer")
show_figure(v.fig, "shoebox_interactive")

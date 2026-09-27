# co2lab-spillpoint's trapsHuginWest: the Hugin West formation of the CO2
# Storage Atlas, where the node-based and the cell-based methods differ most
# in trap volume. The whole formation with the node-based method, then its
# southern part with both methods side by side.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/trapping/traps_hugin_west.jl
#
# Plotting backend: GLMakie (a window per figure, 3D views rotatable with the
# mouse; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, Printf

# --------------------------- 1. The whole formation ---------------------------

# The node-based method (co2lab's trapAnalysis(Gt, false)) follows spill
# paths along the nodes and edges of the caprock surface; the result is then
# given per column.
hugin = atlas_mesh("Huginfmwest")
s = trap_surface(hugin.mesh, hugin.columns)
res = trap_analysis(s; method = :node)
@printf("Hugin West: %d columns, %d traps (node-based)\n", length(s.zc),
        maximum(res.traps))
fig = Figure(size = (1100, 900))
ax = ve_axis3(fig[1, 1], s; view = (-70, 55), exaggeration = 20,
              title = "Hugin West: traps and spill regions (node-based)")
plot_traps!(ax, s, res)
show_figure(fig, "hugin_west_traps")

# ------------ 2. The southern part: node- against cell-based traps ------------

# The cells east of 435 km and south of 6450 km (co2lab's extractSubgrid).
x = tpfv_geometry(hugin.mesh).cell_centroids
south = grid_section_subset(hugin.grid, hugin.mesh,
                            (x[1, :] .> 4.35e5) .& (x[2, :] .< 6.45e6))
s = trap_surface(south.mesh, south.columns)
node = trap_analysis(s; method = :node)
cell = trap_analysis(s; method = :cell)

# The two methods number their traps independently: each node-based trap is
# listed with the cell-based trap sharing most of its columns.
match(t) = argmax([count((node.traps .== t) .& (cell.traps .== u))
                   for u in 1:maximum(cell.traps)])
vn, vc = trap_volumes(s, node), trap_volumes(s, cell)
println("\nnode-based trap: columns, volume (m3) | cell-based trap: columns, volume")
for t in 1:maximum(node.traps)
    u = match(t)
    @printf("%4d: %5d  %9.2e | %4d: %5d  %9.2e\n", t, count(==(t), node.traps),
            vn[t], u, count(==(u), cell.traps), vc[u])
end
@printf("total                %9.2e |                    %9.2e\n", sum(vn), sum(vc))

fig = Figure(size = (1500, 750))
for (k, (name, r)) in enumerate(("node-based" => node, "cell-based" => cell))
    a = ve_axis3(fig[1, k], s; view = (-20, 40), exaggeration = 20,
                 title = "southern Hugin West: $name traps")
    plot_traps!(a, s, r; labels = true)
end
show_figure(fig, "hugin_west_south")

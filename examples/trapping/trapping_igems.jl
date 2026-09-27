# co2lab-spillpoint's trappingIGEMS: capacity estimates and spill-point
# analysis on a surface of the IGEMS data set, which studied how top-surface
# morphology influences CO2 storage capacity. The surface used, OSSNP1 no. 1,
# has offshore sand ridges under marine shale (amplitude up to 20 m, 2-4 km
# wide, 10-60 km long) and no faults; 180 000 columns of 100 m.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/trapping/trapping_igems.jl
#
# Needs the IGEMS surfaces (files.nr.no/igems/surfaces.zip) or the IGEMS
# sample in atlas_data_directory(). Plotting backend: GLMakie (a window per
# figure; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead. co2lab builds the surface transposed
# (x and y swapped); here it keeps the file's orientation.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Printf

# ------------------------ 1. The surface and its traps ------------------------

# The IRAP surface as a one-layer grid, 100 m thick.
grid = igems_mesh("OSSNP1", 1)
s = trap_surface(grid.mesh, grid.columns)

# Traps, spill depths, connections, rivers and spill regions with the
# node-based method (co2lab's trapAnalysis(Gt, false)); the results are
# given per column. A spill region holds every column from which CO2 would
# eventually flow into its trap.
res = trap_analysis(s; method = :node)
ntrap = maximum(res.traps)
@printf("%d columns; %d traps holding %d columns\n", length(s.zc), ntrap,
        count(>(0), res.traps))

# The traps: long narrow ones along the crests of the ridges, and many small
# scattered pockets; then the rivers leaving them.
fig = plot_ve_map(s; regions = res.traps .> 0, region_color = :red,
                  depth_contours = 20, title = "OSSNP1: traps (red) and caprock depth")
show_figure(fig, "igems_traps")
fig = plot_traps(s, res; regions = false, title = "OSSNP1: traps and rivers")
show_figure(fig, "igems_rivers")

# ---------------------------- 2. Trapping capacity ----------------------------

# Pore volume in each trap (porosity 0.2), and the ten largest traps.
volumes = trap_volumes(s, res; poro = fill(0.2, length(s.zc)))
total = sum(volumes)
@printf("\ntotal trapping capacity %.2e m3\n\n", total)
order = sortperm(volumes; rev = true)
println("trap   | volume (m3) | columns")
for t in order[1:10]
    @printf("%6d | %11.3e | %7d\n", t, volumes[t], count(==(t), res.traps))
end
@printf("\nthe ten largest traps hold %.2e m3, %.1f %% of the capacity\n",
        sum(volumes[order[1:10]]), 100 * sum(volumes[order[1:10]]) / total)

# The ten largest traps; all traps coloured by volume; the spill regions.
largest = res.traps .∈ Ref(order[1:10])
fig = Figure(size = (1800, 800))
ax = Axis(fig[1, 1]; title = "the ten largest traps (red)", aspect = DataAspect())
plot_ve_map!(ax, s; regions = res.traps .> 0, region_color = :gray40, depth_contours = 0)
plot_ve_map!(ax, s; regions = largest, region_color = :red, depth_contours = 0)
ax = Axis(fig[1, 2]; title = "trap volume (m³)", aspect = DataAspect())
plot_traps!(ax, s, res; regions = false, rivers = false, trap_values = volumes)
Colorbar(fig[2, 2]; limits = extrema(volumes), vertical = false)
ax = Axis(fig[1, 3]; title = "spill regions", aspect = DataAspect())
plot_traps!(ax, s, res; rivers = false)
show_figure(fig, "igems_capacity")

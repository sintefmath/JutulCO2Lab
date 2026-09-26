# VE columns and the top surface (co2lab's showTopSurfaceGrids.m): how a 3D
# grid becomes a VE grid in JutulCO2Lab, on a 3 x 3 x 3 grid with a column cut
# out and a horizontal hole. co2lab builds a separate 2D top-surface grid
# (Gt) for this; JutulCO2Lab keeps the 3D mesh and groups its cells into
# columns (VEColumns), from which the VE domain and the top surface follow.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/atlas/show_top_surface_grids.jl
#
# Plotting backend: GLMakie (a window per figure; close it to continue) by
# default; JUTULCO2LAB_BACKEND=CairoMakie writes the figures to files instead.

backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using JutulCO2Lab, Jutul, JutulDarcy, Printf

# --------------------------- 1. A small 3D grid -------------------------------

# 3 x 3 x 3 cells of 1 x 1 x 0.5 m from 2 m depth, the depth perturbed by
# 0.125 sin(x y), as a corner-point grid section (depth positive down).
nx, ny, nz = 3, 3, 3
depth(x, y, k) = 2.0 + 0.5 * (k - 1) + 0.125 * sin(x * y)

coord = Float64[]
for j in 0:ny, i in 0:nx
    append!(coord, [i, j, 0.0, i, j, 10.0])           # vertical pillars
end
# each cell's eight corners: corner a of 2n along x lies on pillar 1 + a ÷ 2
zcorn = [depth(a ÷ 2, b ÷ 2, (c + 1) ÷ 2 + (isodd(c) ? 0 : 1))
         for c in 1:2nz for b in 1:2ny for a in 1:2nx]

# Cut out column (1, 1) and, as a horizontal hole, cells (3, 1:2, 2).
cut = fill(false, nx, ny, nz)
cut[1, 1, :] .= true
cut[3, 1:2, 2] .= true
section(active) = Dict("cartDims" => (nx, ny, nz), "COORD" => coord,
                       "ZCORN" => zcorn, "ACTNUM" => active)
full = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(section(fill(true, nx, ny, nz)))
mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(section(.!cut))
println("3D grid: ", number_of_cells(mesh), " of 27 cells active")

fig = Figure(size = (1200, 550))
ax1 = Axis3(fig[1, 1]; title = "cells to cut out (red)", zreversed = true)
plot_cell_data!(ax1, full, Float64[cut[cell_ijk(full, c)...] for c in 1:27];
                colormap = [:lightgray, :red])
ax2 = Axis3(fig[1, 2]; title = "the grid", zreversed = true)
plot_mesh!(ax2, mesh; color = :lightgray)
Jutul.plot_mesh_edges!(ax2, mesh; color = :black)
show_figure(fig, "tsg_grid")

# ---------------------------- 2. Columns ---------------------------------------

# A column is a vertically connected stack of cells in one pillar. The hole in
# pillars (3, 1) and (3, 2) disconnects the cell below it: by default, as in
# co2lab, cells below a hole are discarded (unlike co2lab, a hole is a missing
# face between cells, not a gap in the logical k index).
columns, discarded = ve_columns(mesh)
println("\ncolumns: ", number_of_columns(columns), "; discarded cells: ", discarded)
for c in 1:number_of_columns(columns)
    cells = column_cells(columns, c)
    @printf("  column %d: pillar %s, cells %s (top to bottom)\n", c,
            columns.pillar[c], join(cells, ", "))
end

# The discarded cells are then removed from the mesh (ve_trim) before the
# columns are formed again.
mesh = ve_trim(mesh, discarded)
columns, _ = ve_columns(mesh)
println("after ve_trim: ", number_of_cells(mesh), " cells in ",
        number_of_columns(columns), " columns")

# Alternatively the part below a hole becomes a column of its own, in the same
# pillar ("stacked"); such domains need the hybrid model of roadmap Stage 12.
untrimmed = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(section(.!cut))
stacked, _ = ve_columns(untrimmed; below_holes = :stack)
println("with below_holes = :stack: ", number_of_columns(stacked), " columns")

# ------------------------- 3. The top surface ----------------------------------

# The caprock seen from above: each column's top face. The surface is drawn
# slightly above the grid (0.15 m).
surface = ve_top_surface(mesh, columns)
fig = Figure(size = (900, 750))
ax = Axis3(fig[1, 1]; title = "top surface (blue) over the 3D grid", zreversed = true,
           azimuth = 0.5, elevation = 0.5)
plot_mesh!(ax, mesh; color = (:lightgray, 0.5), transparency = true)
Jutul.plot_mesh_edges!(ax, mesh; color = :gray)
plot_ve_surface!(ax, surface; offset = -0.15)
# each column's number and pillar at its centre
for c in 1:number_of_columns(columns)
    z = sum(surface.z[surface.polygons[c]]) / length(surface.polygons[c]) - 0.3
    text!(ax, Point3f(surface.xc[c], surface.yc[c], z);
          text = "$c: $(columns.pillar[c])", align = (:center, :center))
end
show_figure(fig, "tsg_top_surface")

# ----------------------------- 4. Geometry ------------------------------------

# Per column: caprock depth z_top, height H; per cell (top to bottom): its
# thickness dz and its bottom's depth below the caprock, z.
geo = ve_column_geometry(mesh, columns)
println("\ncolumn  z_top     H     cells: dz / z")
for c in 1:number_of_columns(columns)
    r = columns.pos[c]:(columns.pos[c+1] - 1)
    cells = join((@sprintf("%.2f / %.2f", geo.dz[m], geo.z[m]) for m in r), ",  ")
    @printf("  %d     %.3f  %.3f   %s\n", c, geo.z_top[c], geo.H[c], cells)
end

# The VE domain: one cell per column, rock averaged over the column, and the
# gravity term taken between caprock depths.
ve = ve_domain(mesh, columns; perm_x = fill(1e-13, number_of_cells(mesh)),
               poro = fill(0.2, number_of_cells(mesh)))
println("\nVE domain: ", number_of_cells(physical_representation(ve.domain)),
        " cells, ", size(ve.domain[:neighbors], 2), " interior faces")

fig = Figure(size = (900, 750))
ax = Axis3(fig[1, 1]; title = "caprock (red) and bottom (blue) of each column",
           zreversed = true, azimuth = 0.5, elevation = 0.4)
plot_ve_surface!(ax, surface; color = (:red, 0.3))
scatter!(ax, surface.xc, surface.yc, geo.z_top; color = :red, markersize = 12)
scatter!(ax, surface.xc, surface.yc, geo.z_top .+ geo.H; color = :blue,
         markersize = 12)
show_figure(fig, "tsg_geometry")

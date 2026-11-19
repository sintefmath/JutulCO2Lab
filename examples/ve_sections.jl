# Cross-sections of a VE grid: which cells form VE columns and which are fine
# 3D cells, where pillars are cut into stacked columns, and where the barriers
# lie. A sloping aquifer gets the features that shape its columns: a shale
# lens without net rock (removed: a hole), a tight layer (low vertical
# permeability: fine 3D cells), a leaky shale drape and a sealing layer
# (transmissibility multipliers), and fine 3D pillars around a well site.
# Then the Johansen formation's VE columns, if its data set is present.
#
#   julia --project=examples examples/ve_sections.jl
#
# Self-contained (the grid is built here); no simulation.
#
# Plotting backend: GLMakie (a window per figure; resizing it changes the
# vertical exaggeration) by default; JUTULCO2LAB_BACKEND=CairoMakie writes the
# figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, JutulDarcy

mD = si_unit(:milli) * si_unit(:darcy)

# ------------------------------- 1. The aquifer -------------------------------

# co2lab's sloping aquifer, 3 km long, 500 m wide, 30 m thick: 60 × 5 pillars
# of 20 layers (1.5 m)
grid = sloping_aquifer_grid(; nx = 60, ny = 5, nz = 20, length = 3000.0,
                            width = 500.0, thickness = 30.0)
gross = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
ijk = [cell_ijk(gross, c) for c in 1:number_of_cells(gross)]

# ---------------------------------- 2. Rock -----------------------------------

# Net-to-gross 0.6 at the top to 1 at the bottom, but 0 in a shale lens
# (layers 8-9 under pillars i = 15-25); lateral permeability 200 mD, vertical
# 20 mD, but 0.001 mD in a tight layer (layer 14 under i = 30-40)
lens(i, k) = 15 <= i <= 25 && k in 8:9
ntg = [lens(i, k) ? 0.0 : 0.6 + 0.4 * (k - 1) / 19 for (i, _, k) in ijk]
perm_z = [30 <= i <= 40 && k == 14 ? 1e-3mD : 20mD for (i, _, k) in ijk]
nc = length(ijk)

# The net rock: porosity and lateral permeability times net-to-gross, the
# lens's cells removed from the mesh; the given vertical permeability as it is
net = net_rock(gross; poro = fill(0.25, nc), perm_x = fill(200mD, nc), perm_z, ntg)
mesh = net.mesh
println(nc - number_of_cells(mesh), " cells without net rock removed")

# -------------------------------- 3. Barriers ---------------------------------

# Transmissibility multipliers on the faces between layers, under some
# pillars: a leaky shale drape (0.001, between layers 4 and 5, i = 45-55) and a
# sealing layer (0, between layers 17 and 18, i = 5-12); 1 elsewhere
ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
N = mesh.faces.neighbors
function multiplier(f)
    a, b = ijk[N[f][1]], ijk[N[f][2]]
    a[1:2] == b[1:2] || return 1.0
    k = minmax(a[3], b[3])
    45 <= a[1] <= 55 && k == (4, 5) && return 0.001
    5 <= a[1] <= 12 && k == (17, 18) && return 0.0
    return 1.0
end
mult = multiplier.(1:number_of_faces(mesh))

# The columns: cells, and faces with their multiplier, of vertical
# permeability below 0.1 mD are barriers (threshold chosen here, no default;
# the drape: 0.001 · 20 mD); each stretch between barriers and holes a VE
# column of its own (:stack); the pillars i = 56-58 around a well site fine
# 3D cells
columns, _ = ve_columns(mesh; below_holes = :stack, transmissibility_multiplier = mult,
                        perm_z = net.perm_z, barrier_perm_z = 0.1mD,
                        is_ve = ve_except_pillars(mesh, (i, j) -> 56 <= i <= 58))
println(count(columns.is_ve), " VE columns and ", count(!, columns.is_ve),
        " fine 3D cells in ", length(unique(columns.column_ij)), " pillars")

# ------------------------------ 4. The sections -------------------------------

# Along the aquifer (i) through its middle row j = 3: the lens's hole, the
# tight layer's 3D cells, the drape (dashed) and the sealing layer (solid)
# cutting their pillars, the fine pillars at the well site
along_i = ve_section(mesh, columns; along = :i, index = 3)
show_figure(plot_ve_section(along_i; exaggeration = 10,
                            title = "Along the aquifer, j = 3"), "ve_section_i")

# The same without the fine cells inside the VE columns, barriers in purple
show_figure(plot_ve_section(along_i; exaggeration = 10, internal = nothing,
                            sealed = (color = :purple, linewidth = 3),
                            leaky = (color = :purple, linewidth = 3,
                                     linestyle = :dot),
                            title = "VE columns only, j = 3"), "ve_section_i_blocks")

# Across the aquifer (j) through the lens, at i = 20
across = ve_section(mesh, columns; along = :j, index = 20)
show_figure(plot_ve_section(across; exaggeration = 5, title = "Across, i = 20"),
            "ve_section_j")

# -------------------------------- 5. Johansen ---------------------------------

# The Johansen formation's VE columns (co2lab's VE grid: below the first hole
# discarded, faults between pillars), along i through the injector's row
if isdir(johansen_data_directory())
    johansen = johansen_ve_domain()
    s = ve_section(johansen.mesh, johansen.columns; along = :i, index = 48)
    show_figure(plot_ve_section(s; exaggeration = 5, title = "Johansen, j = 48"),
                "ve_section_johansen")
else
    println("Johansen data set not found: its section skipped")
end

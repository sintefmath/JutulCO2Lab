# The Norwegian CO2 Storage Atlas: a formation as a VE domain (roadmap Stage 5,
# step 3) — co2lab-common/atlasgrid/getFormationTopGrid.m.
#
# co2lab: getAtlasGrid, refineGrdecl (vertical layers), processGRDECL, the
# largest connected piece, topSurfaceGrid, averageRock. Here: the grid section
# of atlas_grid.jl, GeoEnergyIO's mesh_from_grid_section, the largest piece,
# then ve_columns and ve_domain as for any other 3D grid (roadmap Stage 1).

# ------------------------------- The VE domain --------------------------------

"""
    atlas_ve_domain(name; coarsening = 1, layers = 1, min_thickness = 1e-3,
                    perm = nothing, poro = nothing, rugosity = 0.0,
                    temperature = nothing, co2lab_compatible = false) -> NamedTuple

Atlas formation `name` as a VE domain, as co2lab's
`getFormationTopGrid(name, coarsening, 'vertical_layers', layers)`:

1. the grid section of [`atlas_grid_section`](@ref) with `nz = layers` (co2lab
   refines a one-layer deck with `refineGrdecl` instead: the same layers, up
   to rounding);
2. the 3D mesh from GeoEnergyIO's `mesh_from_grid_section`, restricted to its
   largest connected piece (the other pieces are deactivated in `ACTNUM` and
   the mesh built again);
3. columns ([`ve_columns`](@ref)) and the VE domain ([`ve_domain`](@ref));
4. rock as co2lab: with rock maps, the maps per 3D cell, their `NaN`s
   replaced by the formation's averages, then averaged over each column;
   otherwise the formation's average values (`petro`) everywhere. With rock
   maps, net-to-gross comes only from a net-to-gross map (co2lab drops the
   Atlas's average net-to-gross then, and so does this); without, from the
   average if the Atlas gives one.

`perm` (m²) and `poro`, if given, replace the formation's average values:
needed where the Atlas has none (the Brent Group, Pliocene sand, Sleipner,
and the sealing Rør and Not formations), for which co2lab builds grids with
`NaN` rock; here they throw an `ArgumentError` instead. `rugosity` (m, a
number or one value per column) and `temperature` (K; a number, one value
per column or a function of depth) are passed on to [`ve_domain`](@ref);
`min_thickness` (per layer) and `co2lab_compatible` to
[`atlas_grid_section`](@ref).

Returns the NamedTuple `(mesh, columns, ve, grid, petro)`, all new: the 3D
`UnstructuredMesh`, its `VEColumns`, [`ve_domain`](@ref)'s result, the grid
section with `ACTNUM` restricted to the mesh's cells, and the formation's
average rock (`(perm, poro, ntg)`, 40 in `STRUCTURES.md`, with `perm`,
`poro` as given).
"""
function atlas_ve_domain(name::AbstractString; coarsening::Integer = 1,
                         layers::Integer = 1, min_thickness::Real = 1e-3,
                         perm = nothing, poro = nothing, rugosity = 0.0,
                         temperature = nothing, co2lab_compatible::Bool = false)
    m = atlas_mesh(name; coarsening, layers, min_thickness, co2lab_compatible)
    petro = m.petro
    # given values replace the Atlas's averages
    isnothing(perm) || (petro = (; petro..., perm = Float64(perm)))
    isnothing(poro) || (petro = (; petro..., poro = Float64(poro)))
    (isnan(petro.perm) || isnan(petro.poro)) &&
        throw(ArgumentError("the Atlas gives no average permeability or porosity " *
                            "for $name: pass `perm` and `poro`"))

    rock = _atlas_cell_rock(m.grid, m.mesh, petro)
    ve = ve_domain(m.mesh, m.columns; rock..., rugosity, temperature)
    return (mesh = m.mesh, columns = m.columns, ve = ve, grid = m.grid, petro = petro)
end

"""
    atlas_mesh(name; coarsening = 1, layers = 1, min_thickness = 1e-3,
               co2lab_compatible = false) -> NamedTuple

Steps 1–3 of [`atlas_ve_domain`](@ref) without the rock: the geometry of
Atlas formation `name` — also of formations without Atlas rock values (for
maps, trap analysis). The arguments are as for [`atlas_ve_domain`](@ref).

Returns the NamedTuple `(mesh, columns, grid, petro)`, all new: the 3D
`UnstructuredMesh` of the largest connected piece, its `VEColumns`, the grid
section with `ACTNUM` restricted to the mesh's cells, and the Atlas's average
rock (possibly `NaN`).
"""
function atlas_mesh(name::AbstractString; coarsening::Integer = 1,
                    layers::Integer = 1, min_thickness::Real = 1e-3,
                    co2lab_compatible::Bool = false)
    s = atlas_grid_section(name; coarsening, nz = layers, min_thickness,
                           co2lab_compatible)
    grid = s.grid
    # the largest connected piece only, as co2lab's G(1)
    grid["ACTNUM"] = _largest_piece(grid)
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    columns, _ = ve_columns(mesh)
    return (mesh = mesh, columns = columns, grid = grid, petro = s.petro)
end

"""
    grid_section_subset(grid, mesh, keep) -> NamedTuple

Part of a corner-point grid (co2lab's `extractSubgrid`): the grid section
`grid` with only the cells of its mesh `mesh` for which `keep` (one `Bool`
per cell of `mesh`) is true left active, meshed anew — e.g. the cells whose
centroids lie in a region. `grid` and `mesh` are only read.

Returns `(mesh, columns, grid)`, all new: the mesh, its `VEColumns` and the
grid section (a copy of `grid` with a new `"ACTNUM"`).
"""
function grid_section_subset(grid, mesh, keep::AbstractVector{Bool})
    length(keep) == number_of_cells(mesh) ||
        throw(ArgumentError("keep must have one value per cell of mesh"))
    # the kept cells active, all others not; then meshed anew
    active = fill(false, Tuple(grid["cartDims"]))
    for c in findall(keep)
        active[cell_ijk(mesh, c)...] = true
    end
    sub = merge(grid, Dict("ACTNUM" => active))
    m = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(sub)
    return (mesh = m, columns = first(ve_columns(m)), grid = sub)
end

# ACTNUM (new Array{Bool, 3}) keeping only the largest connected piece of the
# section's mesh (co2lab: processGRDECL's pieces, the one with most cells)
function _largest_piece(grid)
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    # pieces connected through faces
    N = get_neighborship(mesh)
    g = Graphs.SimpleGraph(number_of_cells(mesh))
    foreach(f -> Graphs.add_edge!(g, N[1, f], N[2, f]), axes(N, 2))
    piece = argmax(length, Graphs.connected_components(g))

    active = fill(false, grid["cartDims"])
    for c in piece
        active[cell_ijk(mesh, c)...] = true
    end
    return active
end

# ------------------------------------ Rock ------------------------------------

# Per 3D cell of `mesh` (NamedTuple for ve_domain: perm_x, poro, ntg): the
# section's rock maps with NaN replaced by the averages (co2lab's updateNaNs),
# or the averages where there are no maps
function _atlas_cell_rock(grid, mesh, petro)
    nc = number_of_cells(mesh)
    # a map's value per cell, its gaps (and missing maps) filled with the average
    at(key) = [grid[key][cell_ijk(mesh, c)...] for c in 1:nc]
    from_map(key, avg) = haskey(grid, key) ? replace(at(key), NaN => avg) :
                         fill(avg, nc)

    perm = from_map("PERMX", petro.perm)
    poro = from_map("PORO", petro.poro)
    # co2lab: with any rock map, net-to-gross only from a map
    maps = any(k -> haskey(grid, k), ("PERMX", "PORO", "NTG"))
    has_ntg = haskey(grid, "NTG") || (!maps && !isnothing(petro.ntg))
    ntg = has_ntg ? from_map("NTG", something(petro.ntg, NaN)) : nothing
    return (perm_x = perm, poro = poro, ntg = ntg)
end

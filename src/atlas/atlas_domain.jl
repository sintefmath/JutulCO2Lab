# The Norwegian CO2 Storage Atlas, part 3: a formation as a 3D mesh and a VE
# domain (co2lab's getFormationTopGrid).
#
# The grid section of atlas_grid.jl is meshed by GeoEnergyIO, cut to its
# largest connected piece, then given columns and a VE domain (ve_columns,
# ve_domain) as any other 3D grid. Also: subsets of a grid section.

# ------------------------------- The VE domain --------------------------------

"""
    atlas_ve_domain(name; coarsening = 1, layers = 1, min_thickness = 1e-3,
                    perm = nothing, poro = nothing, rugosity = 0.0,
                    temperature = nothing, co2lab_compatible = false) -> NamedTuple

An Atlas formation as a VE domain, ready for
[`setup_co2lab_model`](@ref); co2lab's
`getFormationTopGrid(name, coarsening, 'vertical_layers', layers)`.

Steps: the mesh of [`atlas_mesh`](@ref) (grid section, largest connected
piece, columns), then the rock per 3D cell, its net rock ([`net_rock`](@ref))
and [`ve_domain`](@ref).

# Arguments

- `name`: the formation (case-insensitive).
- `coarsening`, `min_thickness` (m, per layer), `co2lab_compatible`: as
  [`atlas_grid_section`](@ref).
- `layers`: the number of vertical layers of the 3D grid (`nz`). Default 1.
- `perm` (m²), `poro`: replace the formation's Atlas average values. Needed
  where the Atlas has none (the Brent Group, Pliocene sand, Sleipner, and
  the sealing Rør and Not formations); without them those throw an
  `ArgumentError` (co2lab builds grids with `NaN` rock instead).
- `rugosity` (m), `temperature` (K): passed on to [`ve_domain`](@ref); each
  a number or one value per column, `temperature` also a function of depth.

# Returns

The NamedTuple `(mesh, columns, ve, grid, petro, rock)`, all new:

- `mesh`: the 3D `UnstructuredMesh` of the largest connected piece.
- `columns`: its `VEColumns`.
- `ve`: [`ve_domain`](@ref)'s result, on the net rock.
- `grid`: the grid section, `"ACTNUM"` restricted to the mesh's cells.
- `petro`: the formation's average rock `(perm, poro, ntg)`, with `perm` and
  `poro` as given.
- `rock`: the gross rock per cell of `mesh`, `(perm_x, poro, ntg)`, as the
  Atlas gives it (`ntg` is `nothing` without net-to-gross); its
  [`net_rock`](@ref) is what `ve` holds, and what a vertical model over the
  cells (e.g. [`LayeredFringe`](@ref)) takes.

# Notes

- Rock, as co2lab: with rock maps, the maps per 3D cell with their `NaN`s
  replaced by the averages, then averaged over each column; otherwise the
  averages everywhere.
- With rock maps, net-to-gross comes only from a net-to-gross map (co2lab
  then drops the Atlas's average net-to-gross, and so does this); without
  maps, from the average if the Atlas gives one.
- Net-to-gross is folded into porosity and permeability per cell
  ([`net_rock`](@ref)), so `ve.rock` and the domain hold net values where
  co2lab's `rock2D` holds gross ones and a separate `ntg`. Net-to-gross is
  constant over each Atlas column, so pore volumes and transmissibilities are
  co2lab's. Cells with net-to-gross 0 would be removed and the mesh built
  anew; no Atlas map has any.
- co2lab makes the layers by refining a one-layer deck (`refineGrdecl`);
  the layers are the same up to rounding.
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
    # cells without net rock removed, the rest meshed anew
    if !isnothing(rock.ntg) && any(<=(0), rock.ntg)
        m = (; m..., grid_section_subset(m.grid, m.mesh, rock.ntg .> 0)...)
        rock = _atlas_cell_rock(m.grid, m.mesh, petro)
    end

    # the domain on the net rock: net-to-gross in porosity and permeability
    net = net_rock(m.mesh; rock...)
    ve = ve_domain(m.mesh, m.columns; perm_x = net.perm_x, poro = net.poro, rugosity,
                   temperature)
    return (mesh = m.mesh, columns = m.columns, ve = ve, grid = m.grid, petro = petro,
            rock = rock)
end

"""
    atlas_mesh(name; coarsening = 1, layers = 1, min_thickness = 1e-3,
               co2lab_compatible = false) -> NamedTuple

The geometry of an Atlas formation, without rock: the 3D mesh and its
columns, for maps and trap analysis. Works also for formations without Atlas
rock values. [`atlas_ve_domain`](@ref) builds on it.

The grid section of [`atlas_grid_section`](@ref) (with `nz = layers`) is
meshed by GeoEnergyIO and restricted to its largest connected piece (the
other pieces are made inactive and the mesh built again).

# Arguments

As for [`atlas_ve_domain`](@ref): `name`, `coarsening`, `layers`,
`min_thickness`, `co2lab_compatible`.

# Returns

The NamedTuple `(mesh, columns, grid, petro)`, all new:

- `mesh`: the 3D `UnstructuredMesh` of the largest connected piece.
- `columns`: its `VEColumns` ([`ve_columns`](@ref)).
- `grid`: the grid section, `"ACTNUM"` restricted to the mesh's cells.
- `petro`: the Atlas's average rock `(perm, poro, ntg)`, possibly `NaN`.
"""
function atlas_mesh(name::AbstractString; coarsening::Integer = 1,
                    layers::Integer = 1, min_thickness::Real = 1e-3,
                    co2lab_compatible::Bool = false)
    s = atlas_grid_section(name; coarsening, nz = layers, min_thickness,
                           co2lab_compatible)
    grid = s.grid
    # the largest connected piece only, as co2lab's G(1)
    grid["ACTNUM"] = _largest_piece(grid)
    mesh = _mesh_from_grid_section(grid)
    columns, _ = ve_columns(mesh)
    return (mesh = mesh, columns = columns, grid = grid, petro = s.petro)
end

"""
    grid_section_subset(grid, mesh, keep) -> NamedTuple

Part of a corner-point grid, meshed anew (co2lab's `extractSubgrid`), e.g.
the cells whose centroids lie in a region.

# Arguments

- `grid`: a grid section.
- `mesh`: the mesh of `grid`.
- `keep`: one `Bool` per cell of `mesh`, true for the cells to keep. A wrong
  length throws an `ArgumentError`.

# Returns

The NamedTuple `(mesh, columns, grid)`, all new: the mesh of the kept cells,
its `VEColumns`, and a copy of `grid` with a new `"ACTNUM"`.
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
    m = _mesh_from_grid_section(sub)
    return (mesh = m, columns = first(ve_columns(m)), grid = sub)
end

# A new ACTNUM (Array{Bool, 3}) keeping only the largest connected piece of
# the section's mesh (co2lab: the processGRDECL piece with most cells)
function _largest_piece(grid)
    mesh = _mesh_from_grid_section(grid)
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

# Per 3D cell of `mesh` (NamedTuple for net_rock: perm_x, poro, ntg): the
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

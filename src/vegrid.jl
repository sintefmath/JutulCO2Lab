# VE columns on a 3D parent mesh.
#
# The VE domain is a Jutul CoarseMesh over the 3D grid, one block per column
# (see analysis/ve-meshing/grid-representation-comparison.md). CoarseMesh
# accepts any partition and knows nothing about columns; this file builds
# the column partition and the per-column data VE needs, from the 3D mesh
# alone.
#
# Column identification follows MRST's topSurfaceGrid/identify_column_cells
# (co2lab-common/grid/topSurfaceGrid.m:735-779): purely logical. Cells are
# grouped by (i, j) pillar and ordered by k — no geometric sort.

# ------------------------------ Column partition ------------------------------

"""
    VEColumns

Column partition of a 3D mesh.

- `cells`: fine cells of all columns, concatenated, each column top→bottom
- `pos`: column `c` is `cells[pos[c]:pos[c+1]-1]`
- `partition`: fine cell → column, as `CoarseMesh(mesh, partition)` expects
- `pillar`: logical `(i, j)` of each column

Columns are numbered in increasing pillar order (`i + nx*(j-1)`), which is the
order MRST gives the cells of its top-surface grid.
"""
struct VEColumns
    cells::Vector{Int}
    pos::Vector{Int}
    partition::Vector{Int}
    pillar::Vector{Tuple{Int, Int}}
end

"""
    number_of_columns(columns) -> Int

Number of columns in `columns`.
"""
number_of_columns(c::VEColumns) = length(c.pos) - 1

"""
    column_cells(columns, col) -> view

The fine cells of column `col`, top to bottom, as indices into the mesh the
columns were built on. Returns a **view** into `columns.cells`, not a copy:
writing to it would change `columns`.
"""
column_cells(c::VEColumns, col::Integer) = view(c.cells, c.pos[col]:(c.pos[col+1]-1))

# ------------------------------- Mesh from MRST -------------------------------

"""
    ve_mesh_from_mrst(G::AbstractDict) -> UnstructuredMesh

Build a Jutul mesh from an MRST grid read with `MAT.matread`, keeping MRST's
logical indices (`cartDims`, `cells.indexMap`), which Jutul's own
`UnstructuredMesh(::MRSTWrapMesh)` conversion drops, and the direction of
every boundary face (MRST's face tags, stored as entity tags — see
`vefaces.jl`). Depth is positive down, as in MRST.

Returns a new `UnstructuredMesh`: MRST's cells in MRST's order, its interior
and boundary faces each in MRST's face order, `structure = cartDims`,
`cell_map = cells.indexMap`, `z_is_depth = true`, and a `:direction` tag on
every boundary face. `G` is only read.
"""
function ve_mesh_from_mrst(G::AbstractDict)
    cd = Int.(vec(G["cartDims"]))
    length(cd) == 3 || throw(ArgumentError("expected a 3D grid, got cartDims = $cd"))
    mesh = UnstructuredMesh(Jutul.MRSTWrapMesh(G);
        structure = CartesianIndex(cd[1], cd[2], cd[3]),
        cell_map = Int.(vec(G["cells"]["indexMap"])),
        z_is_depth = true)
    return _set_boundary_directions!(mesh, _mrst_boundary_directions(G, mesh))
end

# ---------------------------- Columns and trimming ----------------------------

"""
    ve_columns(mesh; discard_below_holes = true) -> (columns, discarded)

Identify the VE columns of a 3D `mesh` with logical `(i, j, k)` indices.

Cells are grouped by `(i, j)` and ordered by `k` (top→bottom; `k` is assumed
to increase downwards, as in MRST). With `discard_below_holes = true`, every
cell below the first gap in `k` within a column is discarded — MRST's default,
on the grounds that VE is only valid above an impermeable layer.

Returns `(columns, discarded)`: a new `VEColumns` of the cells kept (cell
indices into `mesh`), and the discarded cells as a sorted `Vector{Int}` of
indices into `mesh` (empty if nothing is discarded). `mesh` is only read.

If any cells are discarded, the partition still refers to `mesh`'s cell
numbering and does **not** cover every cell, so it cannot be passed to
`CoarseMesh` as-is; use [`ve_trim`](@ref) first.
"""
function ve_columns(mesh; discard_below_holes::Bool = true)
    nc = number_of_cells(mesh)
    ijk = [cell_ijk(mesh, c) for c in 1:nc]
    nx, ny, _ = grid_dims_ijk(mesh)
    pillars = Dict{Int, Vector{Int}}()
    for c in 1:nc
        i, j, _ = ijk[c]
        push!(get!(pillars, i + nx*(j - 1), Int[]), c)
    end

    cells = Int[]
    pos = [1]
    pillar = Tuple{Int, Int}[]
    discarded = Int[]
    for key in sort!(collect(keys(pillars)))
        col = sort!(pillars[key], by = c -> ijk[c][3])
        nkeep = length(col)
        if discard_below_holes
            for m in 2:length(col)
                if ijk[col[m]][3] != ijk[col[m-1]][3] + 1
                    nkeep = m - 1
                    break
                end
            end
            append!(discarded, col[(nkeep+1):end])
        end
        append!(cells, col[1:nkeep])
        push!(pos, pos[end] + nkeep)
        push!(pillar, ijk[col[1]][1:2])
    end

    partition = zeros(Int, nc)
    for col in 1:(length(pos) - 1), m in pos[col]:(pos[col+1] - 1)
        partition[cells[m]] = col
    end
    return (VEColumns(cells, pos, partition, pillar), sort!(discarded))
end

"""
    ve_trim(mesh, discarded) -> UnstructuredMesh

Remove `discarded` cells from `mesh`, keeping logical indices and boundary
face directions correct. `mesh` is only read.

Returns a new `UnstructuredMesh` with:
- the kept cells, renumbered 1…n in their original order;
- `cell_map` composed with the parent's, so `cell_ijk` still gives each
  cell's original logical index; `structure` and `z_is_depth` unchanged;
- interior faces between two kept cells unchanged; faces between a kept and
  a discarded cell turned into boundary faces; only the nodes still used;
- a `:direction` tag on every boundary face (the new ones from the logical
  indices of the two original cells).

If `discarded` is empty it returns **`mesh` itself**, not a copy.

Works around `Jutul.extract_submesh` setting the submesh's `cell_map` to the
parent's *local* cell numbers rather than composing it with the parent's own
`cell_map` (`Jutul/src/meshes/unstructured/utils.jl`, fixed upstream in
Jutul.jl #357), which gives wrong `(i, j, k)` for any parent that already has
inactive cells.
"""
function ve_trim(mesh::UnstructuredMesh, discarded)
    isempty(discarded) && return mesh
    keep = setdiff(1:number_of_cells(mesh), discarded)
    nc = number_of_cells(mesh)
    parent_map = isnothing(mesh.cell_map) ? collect(1:nc) : mesh.cell_map
    trimmed = extract_submesh(mesh, keep; cell_map = parent_map[keep],
                              z_is_depth = mesh.z_is_depth)
    dirs = _trimmed_boundary_directions(mesh, keep, trimmed)
    return _set_boundary_directions!(trimmed, dirs)
end

# ------------------------------ Column geometry -------------------------------

"""
    ve_column_geometry(mesh, columns) -> NamedTuple

Per-column and per-fine-cell VE geometry, from the 3D mesh alone:

- `z_top`: caprock depth per column — depth of the centroid of the top
  column cell's top face (MRST's `Gt.cells.z`)
- `H`: column height, `sum(dz)` (MRST's `Gt.cells.H`)
- `dz`: per fine cell (in `columns.cells` order), depth difference between the
  centroids of its bottom and top faces (MRST's `Gt.columns.dz`)
- `z`: per fine cell, cumulative `dz` down the column — depth of the cell's
  bottom face relative to the caprock (MRST's `Gt.columns.z`)
- `dx`, `dy`: per column, the x and y extent of the bounding box of the top
  cell's top face nodes (MRST's `cellDims(Gt)`, used for the well index)
- `z_cell`: per fine cell (in `columns.cells` order), depth of the cell
  centroid (MRST's `G.cells.centroids(:, 3)`, used for 3D pressure)

All seven are new `Vector{Float64}`s, returned as the NamedTuple
`(z_top, H, dz, z, dx, dy, z_cell)`.

A cell's top and bottom faces are found by topology, as MRST's face tags 5
and 6: interior faces from the logical indices of the two cells, boundary
faces from the tags stored at import (`vefaces.jl`). Each cell must have
exactly one of each. `mesh` and `columns` are only read.
"""
function ve_column_geometry(mesh, columns::VEColumns)
    geo = tpfv_geometry(mesh)
    nc = number_of_cells(mesh)
    bdir = _boundary_directions(mesh)
    ijk = [cell_ijk(mesh, c) for c in 1:nc]
    ztop_cell = zeros(nc)
    zbot_cell = zeros(nc)
    for c in 1:nc
        ztop_cell[c], zbot_cell[c] = _cell_top_bottom_depth(mesh, geo, ijk, bdir, c)
    end

    ncol = number_of_columns(columns)
    dz = [zbot_cell[c] - ztop_cell[c] for c in columns.cells]
    z = similar(dz)
    H, z_top, dx, dy = zeros(ncol), zeros(ncol), zeros(ncol), zeros(ncol)
    for col in 1:ncol
        r = columns.pos[col]:(columns.pos[col+1] - 1)
        z[r] = cumsum(dz[r])
        H[col] = sum(dz[r])
        top = columns.cells[first(r)]
        z_top[col] = ztop_cell[top]
        dx[col], dy[col] = _top_face_xy_extent(mesh, ijk, bdir, top)
    end

    z_cell = [geo.cell_centroids[3, c] for c in columns.cells]
    return (z_top = z_top, H = H, dz = dz, z = z, dx = dx, dy = dy, z_cell = z_cell)
end

"""
    _cell_top_bottom_depth(mesh, geo, ijk, bdir, c) -> (z_top, z_bottom)

Centroid depths of cell `c`'s top and bottom faces, identified by direction.
"""
function _cell_top_bottom_depth(mesh, geo, ijk, bdir, c)
    top, bot = Float64[], Float64[]
    for f in mesh.faces.cells_to_faces[c]
        d, _ = _interior_face_direction(mesh, ijk, c, f)
        d == _TOP && push!(top, geo.face_centroids[3, f])
        d == _BOTTOM && push!(bot, geo.face_centroids[3, f])
    end

    for f in mesh.boundary_faces.cells_to_faces[c]
        bdir[f] == _TOP && push!(top, geo.boundary_centroids[3, f])
        bdir[f] == _BOTTOM && push!(bot, geo.boundary_centroids[3, f])
    end

    length(top) == 1 && length(bot) == 1 ||
        error("cell $c: expected one top and one bottom face, " *
              "found $(length(top)) and $(length(bot))")
    return (top[1], bot[1])
end

"""
    _top_face_xy_extent(mesh, ijk, bdir, c) -> (dx, dy)

x and y extent of the bounding box of the nodes of cell `c`'s top face.
"""
function _top_face_xy_extent(mesh, ijk, bdir, c)
    for f in mesh.boundary_faces.cells_to_faces[c]
        bdir[f] == _TOP && return _xy_extent(mesh, mesh.boundary_faces.faces_to_nodes[f])
    end

    for f in mesh.faces.cells_to_faces[c]
        d, _ = _interior_face_direction(mesh, ijk, c, f)
        d == _TOP && return _xy_extent(mesh, mesh.faces.faces_to_nodes[f])
    end
    error("cell $c: no top face")
end

# x and y extent of the bounding box of the given nodes
function _xy_extent(mesh, nodes)
    x = [mesh.node_points[n][1] for n in nodes]
    y = [mesh.node_points[n][2] for n in nodes]
    return (maximum(x) - minimum(x), maximum(y) - minimum(y))
end

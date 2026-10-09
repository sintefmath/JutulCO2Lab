# VE columns on a 3D parent mesh.
#
# The VE domain is a Jutul CoarseMesh over the 3D grid, one block per column
# (see analysis/ve-meshing/grid-representation-comparison.md). CoarseMesh
# accepts any partition and knows nothing about columns; this file builds
# the column partition and the per-column data VE needs, from the 3D mesh
# alone.
#
# A column is a vertically connected stretch of one pillar — in pure VE
# usually the whole pillar (analysis/ve-meshing/isve-branch-proposal.md).
# Cells are grouped by (i, j) pillar and ordered by k, as MRST's
# topSurfaceGrid/identify_column_cells (co2lab-common/grid/topSurfaceGrid.m:
# 735-779) — no geometric sort. Unlike co2lab, a pillar is cut where two
# consecutive cells share no face (or only a sealed one), not at every gap in
# k: a pinched-out layer does not cut (analysis/co2lab-issues.org).

# ------------------------------ Column partition ------------------------------

"""
    VEColumns

Column partition of a 3D mesh. A column is a vertically connected stretch of
one pillar; with `below_holes = :stack` ([`ve_columns`](@ref)) a pillar can
hold several. On a hybrid grid (roadmap Stage 12) a "column" is a *block*: a
VE block (a vertically connected stretch of VE cells) or a single fine 3D
cell.

- `cells`: fine cells of all columns, concatenated, each column top→bottom
- `pos`: column `c` is `cells[pos[c]:pos[c+1]-1]`
- `partition`: fine cell → column, as `CoarseMesh(mesh, partition)` expects
- `column_ij`: per column, the logical `(i, j)` of its pillar (one entry per
  column, so a pillar holding several columns appears several times)
- `is_ve`: per column, `true` for a VE block, `false` for a fine 3D cell
- `sealed`, `leaky`: the interior faces of the mesh given as sealed (no
  flow) and as leaky (cut into blocks, flow kept), sorted

Columns are numbered in increasing pillar order (`i + nx*(j-1)`), top to bottom
within a pillar — with one column per pillar the order MRST gives the cells
of its top-surface grid.
"""
struct VEColumns
    cells::Vector{Int}                  # fine cells, column by column, top first
    pos::Vector{Int}                    # column c: cells[pos[c]:pos[c+1]-1]
    partition::Vector{Int}              # fine cell → column
    column_ij::Vector{Tuple{Int, Int}}  # per column: its pillar's logical (i, j)
    is_ve::Vector{Bool}                 # per column: VE block (false: one fine cell)
    sealed::Vector{Int}                 # sealed interior faces (no flow)
    leaky::Vector{Int}                  # leaky interior faces (cut, flow kept)
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
    ve_columns(mesh; below_holes = :discard, sealed_faces = (), leaky_faces = (),
               is_ve = nothing) -> (columns, discarded)

Identify the VE columns of a 3D `mesh` with logical `(i, j, k)` indices.

Cells are grouped by `(i, j)` pillar and ordered by `k` (top→bottom; `k` is
assumed to increase downwards, as in MRST). A pillar is cut wherever two
consecutive cells are *disconnected*: they share no interior face, or only
faces listed in `sealed_faces` or `leaky_faces` (interior face indices of
`mesh`; leaky faces separate VE blocks but keep their flow, roadmap Stage
12). A gap in `k` whose cells still share a face (a pinch-out) does not cut.
Below a cut:

- `below_holes = :discard` — co2lab's treatment (`topSurfaceGrid`'s
  `discard_below_holes`): every cell below the pillar's first cut is
  discarded; one column per pillar.
- `below_holes = :stack` — each connected stretch is its own column. Such
  columns connect laterally through partial faces, which take the transition
  flux (roadmap Stage 12, `transition_flux.jl`).

`is_ve` (one `Bool` per cell of `mesh`, `nothing` = all VE) marks fine 3D
cells (`false`) of a hybrid grid (roadmap Stage 12; whole pillars with
[`ve_except_pillars`](@ref)): each becomes a block of its own, and a
pillar's VE cells are also cut where fine cells interrupt them (not a hole:
nothing is discarded for it).

Returns `(columns, discarded)`: a new `VEColumns` of the cells kept (cell
indices into `mesh`), and the discarded cells as a sorted `Vector{Int}` of
indices into `mesh` (empty if nothing is discarded). `mesh` is only read.

If any cells are discarded, the partition still refers to `mesh`'s cell
numbering and does **not** cover every cell, so it cannot be passed to
`CoarseMesh` as-is; use [`ve_trim`](@ref) first.
"""
function ve_columns(mesh; below_holes::Symbol = :discard, sealed_faces = (),
                    leaky_faces = (), is_ve = nothing)
    below_holes in (:discard, :stack) ||
        throw(ArgumentError("below_holes must be :discard or :stack, got $below_holes"))
    nc = number_of_cells(mesh)
    ve = _cell_is_ve(is_ve, nc)
    sealed, leaky = sort!(unique(collect(Int, sealed_faces))),
                    sort!(unique(collect(Int, leaky_faces)))
    ijk = [cell_ijk(mesh, c) for c in 1:nc]
    # per cell, the cell it connects to directly below (0: none)
    below = _connected_below(mesh, ijk, union(sealed, leaky))

    cells, pos, pillar, block_ve, discarded = Int[], [1], Tuple{Int, Int}[], Bool[],
                                              Int[]
    for col in _pillar_cells(mesh, ijk)
        # a pillar falls apart where a cell is not connected to the one below
        stretches = _connected_stretches(col, below)
        # :discard keeps each pillar's top stretch only (co2lab's default)
        if below_holes == :discard
            foreach(r -> append!(discarded, col[r]), stretches[2:end])
            stretches = stretches[1:1]
        end

        # each stretch in VE blocks, cut where fine cells interrupt it
        for r in stretches, (b, isve) in _split_at_fine(col, r, ve)
            append!(cells, col[b])
            push!(pos, pos[end] + length(b))
            push!(pillar, ijk[col[1]][1:2])
            push!(block_ve, isve)
        end
    end

    # fine cell → column (0 for discarded cells)
    partition = zeros(Int, nc)
    for c in 1:(length(pos) - 1), m in pos[c]:(pos[c+1] - 1)
        partition[cells[m]] = c
    end
    columns = VEColumns(cells, pos, partition, pillar, block_ve, sealed, leaky)
    return (columns, sort!(discarded))
end

# One Bool per cell: true for VE cells (all for `nothing`)
function _cell_is_ve(is_ve, nc)
    isnothing(is_ve) && return fill(true, nc)
    length(is_ve) == nc || throw(ArgumentError("is_ve: one value per cell expected"))
    return collect(Bool, is_ve)
end

# The blocks of the stretch `r` of a pillar's cells `col`: runs of VE cells,
# and each fine cell alone; (index range into `col`, is VE) pairs, top first
function _split_at_fine(col, r, ve)
    out = Tuple{UnitRange{Int}, Bool}[]
    first = 0
    for m in r
        if ve[col[m]]
            first == 0 && (first = m)
            continue
        end
        first > 0 && push!(out, (first:(m - 1), true))
        push!(out, (m:m, false))
        first = 0
    end
    first > 0 && push!(out, (first:last(r), true))
    return out
end

"""
    ve_except_pillars(mesh, fine) -> Vector{Bool}

The `is_ve` argument of [`ve_columns`](@ref) for a hybrid grid whose fine
3D regions are whole pillars: `fine` gives the pillars as logical `(i, j)`
pairs (any collection) or as a function `(i, j) -> Bool`.

Returns a new `Vector{Bool}`, one entry per cell of `mesh`: `false` for the
cells of the chosen pillars, `true` elsewhere. `mesh` is only read.
"""
function ve_except_pillars(mesh, fine)
    # a collection of pillars: membership; a function: its value
    isfine = fine isa Function ? ij -> fine(ij...) :
             (set = Set(Tuple{Int, Int}.(fine)); ij -> ij in set)
    return [!isfine(cell_ijk(mesh, c)[1:2]) for c in 1:number_of_cells(mesh)]
end

"""
    _pillar_cells(mesh, ijk) -> Vector{Vector{Int}}

The cells of each pillar that has any, sorted by `k`, pillars in increasing
order `i + nx*(j-1)`.
"""
function _pillar_cells(mesh, ijk)
    nx = grid_dims_ijk(mesh)[1]
    pillars = Dict{Int, Vector{Int}}()
    for c in eachindex(ijk)
        i, j, _ = ijk[c]
        # pillar number with i fastest, MRST's column order
        push!(get!(pillars, i + nx*(j - 1), Int[]), c)
    end
    keys_sorted = sort!(collect(keys(pillars)))
    return [sort!(pillars[key], by = c -> ijk[c][3]) for key in keys_sorted]
end

"""
    _connected_below(mesh, ijk, sealed_faces) -> Vector{Int}

For every cell, the cell below it in its pillar that it shares an unsealed
interior face with (0 if none). Throws if a cell connects to two cells below.
"""
function _connected_below(mesh, ijk, sealed_faces)
    sealed = Set{Int}(sealed_faces)
    below = zeros(Int, length(ijk))
    for f in 1:number_of_faces(mesh)
        f in sealed && continue
        l, r = mesh.faces.neighbors[f]
        # vertical connections only: both cells on the same pillar
        ijk[l][1:2] == ijk[r][1:2] || continue

        up, down = ijk[l][3] < ijk[r][3] ? (l, r) : (r, l)
        below[up] in (0, down) ||
            error("cell $up connects to two cells below it ($(below[up]), $down)")
        below[up] = down
    end
    return below
end

"""
    _connected_stretches(col, below) -> Vector{UnitRange{Int}}

Split the sorted cells `col` of one pillar where consecutive cells are not
connected; returns index ranges into `col`, top to bottom.
"""
function _connected_stretches(col, below)
    out, first = UnitRange{Int}[], 1
    for m in 2:length(col)
        if below[col[m-1]] != col[m]
            push!(out, first:(m - 1))
            first = m
        end
    end
    push!(out, first:length(col))
    return out
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
    # the kept cells keep their logical indices
    parent_map = isnothing(mesh.cell_map) ? collect(1:nc) : mesh.cell_map
    trimmed = extract_submesh(mesh, keep; cell_map = parent_map[keep],
                              z_is_depth = mesh.z_is_depth)
    # extract_submesh drops the direction tags: carry them over
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
    bdir = _boundary_directions(mesh)
    ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
    # each cell's top and bottom depth (face centroids)
    ztop_cell, zbot_cell = _fine_cell_depths(mesh, geo)

    ncol = number_of_columns(columns)
    dz = [zbot_cell[c] - ztop_cell[c] for c in columns.cells]
    z = similar(dz)
    H, z_top, dx, dy = zeros(ncol), zeros(ncol), zeros(ncol), zeros(ncol)
    for col in 1:ncol
        r = columns.pos[col]:(columns.pos[col+1] - 1)
        # depth of each cell's bottom below the caprock
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
    _fine_cell_depths(mesh, geo = tpfv_geometry(mesh)) -> (z_top, z_bottom)

Every cell's top and bottom depth: the centroid depths of its top and bottom
faces. Returns two new vectors, one entry per cell of `mesh`.
"""
function _fine_cell_depths(mesh, geo = tpfv_geometry(mesh))
    nc = number_of_cells(mesh)
    bdir = _boundary_directions(mesh)
    ijk = [cell_ijk(mesh, c) for c in 1:nc]
    zt, zb = zeros(nc), zeros(nc)
    for c in 1:nc
        zt[c], zb[c] = _cell_top_bottom_depth(mesh, geo, ijk, bdir, c)
    end
    return (zt, zb)
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

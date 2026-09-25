# Face directions and column sides — topology only, no geometric thresholds.
#
# Every fine face has a direction as seen from its cell, in Jutul's naming
# (GeoEnergyIO's corner-point processing, entity tag group `:direction` on
# `BoundaryFaces()`): `:left`/`:right` (logical i−/i+), `:lower`/`:upper`
# (j−/j+), `:top`/`:bottom` (k−/k+). Numbered 1..6 in that order, which is
# also MRST's face-tag order (`G.cells.faces(:, 2)`: W, E, S, N, T, B).
#
# - Interior faces: the direction follows from the logical (i, j, k) of the
#   two cells (`_logical_direction`); the package tests check this against
#   MRST's tags for every interior cell-face of Johansen.
# - Boundary faces have only one cell. Their direction is read from the
#   mesh's `:direction` entity tags: Jutul sets them when it builds a
#   corner-point grid itself; `ve_mesh_from_mrst` sets them from MRST's face
#   tags; `ve_trim` carries them over. Connectivity alone cannot provide it:
#   at a fault, pieces of one side of a cell share edges just like two
#   different sides do.
#
# A column side is the upscaled counterpart of a column: all fine faces of
# the column's cells with the same lateral direction. Two sides are
# neighbours if any of their fine faces is shared; each side has at most one
# neighbour (asserted — true for corner-point grids, whose faults run along
# pillars).

"Face directions, Jutul's names; the numbering 1..6 is MRST's face-tag order."
const _DIRECTIONS = (:left, :right, :lower, :upper, :top, :bottom)
const _TOP = Int8(5)
const _BOTTOM = Int8(6)

"""
    _logical_direction(a, b) -> Int8

Direction (1..6, MRST tag order) of the face between cells with logical
indices `a = (i, j, k)` and `b`, as seen from `a`. Same `(i, j)`: top or
bottom by `k` (increasing downwards). Otherwise exactly one of `i`, `j` must
differ, by one; `k` may differ across a fault.
"""
function _logical_direction(a, b)
    di, dj, dk = b[1] - a[1], b[2] - a[2], b[3] - a[3]
    if di == 0 && dj == 0
        dk == 0 && error("cells $a and $b have the same logical index")
        return dk < 0 ? _TOP : _BOTTOM
    end
    di == -1 && dj == 0 && return Int8(1)
    di == 1 && dj == 0 && return Int8(2)
    di == 0 && dj == -1 && return Int8(3)
    di == 0 && dj == 1 && return Int8(4)
    error("cells $a and $b share a face but are not logically adjacent")
end

"""
    _interior_face_direction(mesh, ijk, c, f) -> (Int8, other_cell)

Direction of interior face `f` as seen from its cell `c`, and the cell on the
other side. `ijk` holds the logical index of every cell of `mesh`.
"""
function _interior_face_direction(mesh, ijk, c, f)
    l, r = mesh.faces.neighbors[f]
    o = l == c ? r : l
    return (_logical_direction(ijk[c], ijk[o]), o)
end

"""
    _boundary_directions(mesh) -> Vector{Int8}

Direction of every boundary face, from the mesh's `:direction` entity tags
(set by Jutul's corner-point processing, [`ve_mesh_from_mrst`](@ref) or
[`ve_trim`](@ref)).
"""
function _boundary_directions(mesh)
    dirs = zeros(Int8, length(mesh.boundary_faces.neighbors))
    for (k, name) in enumerate(_DIRECTIONS)
        faces = Jutul.get_mesh_entity_tag(mesh, BoundaryFaces(), :direction, name; throw = false)
        ismissing(faces) || (dirs[faces] .= k)
    end
    any(==(0), dirs) && error("mesh has boundary faces without a :direction tag (Jutul corner-point grids and ve_mesh_from_mrst / ve_trim provide them)")
    return dirs
end

"""
    _set_boundary_directions!(mesh, dirs)

Store one direction per boundary face as entity tags.

Modifies: `mesh` (its entity tags, group `:direction` on `BoundaryFaces()`).
`dirs` is only read.
"""
function _set_boundary_directions!(mesh, dirs)
    length(dirs) == length(mesh.boundary_faces.neighbors) || error("one direction per boundary face expected")
    for (k, name) in enumerate(_DIRECTIONS)
        Jutul.set_mesh_entity_tag!(mesh, BoundaryFaces(), :direction, name, findall(==(k), dirs))
    end
    return mesh
end

"""
    _mrst_boundary_directions(G, mesh) -> Vector{Int8}

MRST's tag of every boundary face of `mesh`, which must be the direct
conversion of the MRST grid `G`. Jutul keeps MRST's face order when it
separates boundary faces (`UnstructuredMesh`'s raw constructor); that is
checked here exactly — same cell, same node set — for every face.
"""
function _mrst_boundary_directions(G, mesh)
    N = Int.(G["faces"]["neighbors"])
    cf, fpos = G["cells"]["faces"], Int.(vec(G["cells"]["facePos"]))
    fnodes, npos = Int.(vec(G["faces"]["nodes"])), Int.(vec(G["faces"]["nodePos"]))
    tag = zeros(Int8, size(N, 1))
    for c in 1:(length(fpos) - 1), p in fpos[c]:(fpos[c+1] - 1)
        f = Int(cf[p, 1])
        (N[f, 1] == 0 || N[f, 2] == 0) && (tag[f] = Int8(cf[p, 2]))
    end
    bfaces = [f for f in 1:size(N, 1) if N[f, 1] == 0 || N[f, 2] == 0]
    length(bfaces) == length(mesh.boundary_faces.neighbors) || error("boundary face count differs from the MRST grid")
    for (k, f) in enumerate(bfaces)
        same = mesh.boundary_faces.neighbors[k] == N[f, 1] + N[f, 2] &&
               sort(collect(mesh.boundary_faces.faces_to_nodes[k])) == sort(fnodes[npos[f]:(npos[f+1] - 1)])
        same || error("boundary face $k of the mesh is not MRST face $f")
    end
    return tag[bfaces]
end

"""
    _trimmed_boundary_directions(old, keep, new) -> Vector{Int8}

Directions of the boundary faces of `new = extract_submesh(old, keep)`.
Jutul lists them as: interior faces of `old` with exactly one kept cell (in
`old`'s face order), then `old`'s boundary faces with a kept cell. The first
kind get their direction from the logical indices of the two cells; the second
keep theirs. Checked exactly: same cell and same node coordinates per face.
"""
function _trimmed_boundary_directions(old, keep, new)
    active = falses(number_of_cells(old)); active[keep] .= true
    renum = zeros(Int, number_of_cells(old)); renum[keep] = 1:length(keep)
    ijk(c) = cell_ijk(old, c)
    dirs, cells, nodes = Int8[], Int[], Vector{Int}[]
    for f in 1:number_of_faces(old)
        l, r = old.faces.neighbors[f]
        active[l] == active[r] && continue
        c, o = active[l] ? (l, r) : (r, l)
        push!(dirs, _logical_direction(ijk(c), ijk(o))); push!(cells, c); push!(nodes, collect(old.faces.faces_to_nodes[f]))
    end
    olddirs = _boundary_directions(old)
    for bf in eachindex(olddirs)
        c = old.boundary_faces.neighbors[bf]
        active[c] || continue
        push!(dirs, olddirs[bf]); push!(cells, c); push!(nodes, collect(old.boundary_faces.faces_to_nodes[bf]))
    end
    length(dirs) == length(new.boundary_faces.neighbors) || error("unexpected boundary face count after trimming")
    for k in eachindex(dirs)
        same = new.boundary_faces.neighbors[k] == renum[cells[k]] &&
               sort(new.node_points[collect(new.boundary_faces.faces_to_nodes[k])]) == sort(old.node_points[nodes[k]])
        same || error("boundary face $k after trimming does not match the expected face")
    end
    return dirs
end

# ------------------------------------------------------------------ column sides

"""
    VEColumnSides

The lateral sides of every column, one entry per side that has any fine face:

- `column`, `direction` (1..4 = left, right, lower, upper: logical i−, i+, j−, j+)
- `neighbour`: the column on the other side, or 0 — an outer boundary side
- `interior_faces`: fine interior faces shared with `neighbour` (together they
  are the `CoarseMesh`'s coarse face between the two columns)
- `boundary_faces`: fine boundary faces on this side — for a side with a
  neighbour, the exposed parts (e.g. above a lower neighbour), which stay
  no-flow
"""
struct VEColumnSides
    column::Vector{Int}
    direction::Vector{Int8}
    neighbour::Vector{Int}
    interior_faces::Vector{Vector{Int}}
    boundary_faces::Vector{Vector{Int}}
end

Base.length(s::VEColumnSides) = length(s.column)

"""
    ve_column_sides(mesh, columns) -> VEColumnSides

Group the lateral fine faces of every column into sides by direction, and
find each side's neighbouring column from the faces' cells. Throws if a side
has fine faces towards more than one other column. `mesh` must carry face
directions ([`ve_mesh_from_mrst`](@ref), [`ve_trim`](@ref)) and every cell
must belong to a column. `mesh` and `columns` are only read.
"""
function ve_column_sides(mesh, columns::VEColumns)
    bdir = _boundary_directions(mesh)
    ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
    out = VEColumnSides(Int[], Int8[], Int[], Vector{Int}[], Vector{Int}[])
    for col in 1:number_of_columns(columns)
        acc = [(nbrs = Int[], int = Int[], bnd = Int[]) for _ in 1:4]
        for c in column_cells(columns, col)
            _collect_side_faces!(acc, mesh, columns, ijk, bdir, c)
        end
        for d in 1:4
            a = acc[d]
            isempty(a.int) && isempty(a.bnd) && continue
            nb = unique(a.nbrs)
            length(nb) <= 1 || error("column $col, $(_DIRECTIONS[d]) side: faces towards columns $nb; " *
                                     "expected at most one neighbouring column")
            push!(out.column, col); push!(out.direction, d); push!(out.neighbour, isempty(nb) ? 0 : nb[1])
            push!(out.interior_faces, a.int); push!(out.boundary_faces, a.bnd)
        end
    end
    return out
end

"""
    _collect_side_faces!(acc, mesh, columns, ijk, bdir, c)

Add cell `c`'s lateral faces to the per-direction accumulators of its column.

Modifies: `acc` (its `nbrs`, `int`, `bnd` lists). The other arguments are only
read.
"""
function _collect_side_faces!(acc, mesh, columns, ijk, bdir, c)
    col = columns.partition[c]
    for f in mesh.faces.cells_to_faces[c]
        d, o = _interior_face_direction(mesh, ijk, c, f)
        oc = columns.partition[o]
        oc == 0 && error("cell $o is in no column; trim the mesh first (ve_trim)")
        if oc == col
            d in (_TOP, _BOTTOM) || error("lateral face $f inside column $col")
            continue
        end
        d in (_TOP, _BOTTOM) && error("vertical face $f between columns $col and $oc")
        push!(acc[d].int, f); push!(acc[d].nbrs, oc)
    end
    for f in mesh.boundary_faces.cells_to_faces[c]
        bdir[f] <= 4 && push!(acc[bdir[f]].bnd, f)
    end
    return acc
end

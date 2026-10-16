# Vertical cross-sections of a VE or hybrid grid: the cells of one logical row
# of the 3D mesh (fixed j) or one logical column (fixed i) as polygons in the
# section's plane, and their faces seen side-on as line segments by kind
# (block outlines, fine cells inside VE blocks, barriers). Geometry only;
# plotting is in the Makie extension.

# ---------------------------------- Sections ----------------------------------

"""
    ve_section(mesh, columns; along = :i, index) -> NamedTuple

A vertical cross-section through a VE or hybrid grid: the cells of one
logical row of the 3D mesh as polygons in the section's plane, and their faces
seen side-on as segments, by kind. The basis for
[`plot_ve_section`](@ref).

# Arguments

- `mesh`, `columns`: the 3D mesh and its column partition
  ([`ve_columns`](@ref)).
- `along`: `:i` or `:j`, the logical direction the section runs in.
- `index`: the other logical index, fixed (`j` for `along = :i`, `i` for
  `along = :j`).

# Returns

A NamedTuple of new arrays:

- `cells`: the section's cells, as indices into `mesh`, in mesh order; only
  cells in a column (block) of `columns`.
- `polygons`: per cell, its corners `(s, z)` (m): top left, top right,
  bottom right, bottom left. `s` is the distance along the section from its
  start (the smallest `s` is 0), `z` the depth (positive down).
- `block`: per cell, its column (block) in `columns`.
- `is_ve`: per cell, whether its block is a VE block (else a fine 3D cell).
- `segments`: the cells' top, bottom and lateral (along the section) faces
  seen side-on, `((s1, z1), (s2, z2))`, each face once.
- `kind`: per segment, `:outline` (a block's boundary: between two blocks,
  or a boundary face), `:internal` (between two cells of one VE block),
  `:sealed` or `:leaky` (barrier faces, see [`VEColumns`](@ref)).
- `direction`, `origin`: the section's horizontal direction, a unit vector
  `(dx, dy)`, and the value `s0` with `s = dx·x + dy·y − s0`.

Throws an `ArgumentError` for an unknown `along` or a section without cells.

# Notes

- A cell is drawn from its top and bottom faces seen side-on: each face's
  corners are split into the half nearer the section's start and the half
  farther along, and each half averaged (the row's front and back). Lateral
  faces are split by depth the same way. Only a face's corners count: nodes
  on its straight edges (where faces of neighbouring columns do not match)
  are left out.
- The direction runs from the row's first to its last pillar (the mean of
  their cells' top-face nodes), so `s` is close to x (or y) on grids aligned
  with the axes; with a single pillar, the x (or y) axis.
- Faces across the section (the row's front and back) are not segments.
"""
function ve_section(mesh, columns::VEColumns; along::Symbol = :i, index::Integer)
    along in (:i, :j) || throw(ArgumentError("along must be :i or :j, got $along"))
    nc = number_of_cells(mesh)
    ijk = [cell_ijk(mesh, c) for c in 1:nc]
    # the row: the other logical index fixed; cells in a block only
    fixed = along == :i ? 2 : 1
    cells = [c for c in 1:nc if ijk[c][fixed] == index && columns.partition[c] > 0]
    isempty(cells) &&
        throw(ArgumentError("no cells in the section at " *
                            "$(along == :i ? :j : :i) = $index"))

    # node positions in the section's plane: along it, and depth
    bdir = _boundary_directions(mesh)
    d = _section_direction(mesh, ijk, bdir, cells, along)
    pts = mesh.node_points
    sz(n) = (d[1] * pts[n][1] + d[2] * pts[n][2], Float64(pts[n][3]))

    polygons = [_cell_polygon(mesh, ijk, bdir, c, sz) for c in cells]
    segments, kind = _section_segments(mesh, ijk, bdir, columns, cells, along, sz)

    # distances from the section's start
    s0 = minimum(p[1] for poly in polygons for p in poly)
    shift(p) = (p[1] - s0, p[2])
    polygons = [map(shift, poly) for poly in polygons]
    segments = [map(shift, seg) for seg in segments]
    block = columns.partition[cells]
    return (cells = cells, polygons = polygons, block = block,
            is_ve = columns.is_ve[block], segments = segments, kind = kind,
            direction = d, origin = s0)
end

# The section's horizontal unit vector: from the row's first pillar to its last
# (the mean of their cells' top-face nodes); a single pillar: the x or y axis
function _section_direction(mesh, ijk, bdir, cells, along)
    a = along == :i ? 1 : 2
    lo, hi = extrema(ijk[c][a] for c in cells)
    lo == hi && return along == :i ? (1.0, 0.0) : (0.0, 1.0)
    pts = mesh.node_points
    function centre(k)
        ns = [n for c in cells if ijk[c][a] == k
               for n in _face_nodes_toward(mesh, ijk, bdir, c, _TOP)]
        return (sum(n -> pts[n][1], ns), sum(n -> pts[n][2], ns)) ./ length(ns)
    end
    (x0, y0), (x1, y1) = centre(lo), centre(hi)
    L = hypot(x1 - x0, y1 - y0)
    return ((x1 - x0) / L, (y1 - y0) / L)
end

# Cell c's corners in the section's plane: its top and bottom faces side-on,
# (top left, top right, bottom right, bottom left)
function _cell_polygon(mesh, ijk, bdir, c, sz)
    top = _face_corners(mesh, _face_nodes_toward(mesh, ijk, bdir, c, _TOP))
    bottom = _face_corners(mesh, _face_nodes_toward(mesh, ijk, bdir, c, _BOTTOM))
    tl, tr = _split_ends(sz.(top), 1)
    bl, br = _split_ends(sz.(bottom), 1)
    return (tl, tr, br, bl)
end

# A face's corners: its nodes that do not lie between two others. Where
# neighbouring columns' faces do not match, the mesh has extra nodes on a
# face's straight edges (not always listed in order along the edge); averaged
# in, they would pull the face's ends seen side-on towards its middle. All
# nodes if fewer than three corners remain (a degenerate face).
function _face_corners(mesh, nodes)
    pts = mesh.node_points
    # one node per point: a point listed twice would count as between itself
    v = unique(n -> pts[n], nodes)
    between(k) = any(_on_segment(pts[v[a]], pts[v[k]], pts[v[b]])
                     for a in eachindex(v), b in eachindex(v) if a != k != b != a)
    corners = [v[k] for k in eachindex(v) if !between(k)]
    return length(corners) >= 3 ? corners : collect(nodes)
end

# Whether point b lies on the segment from a to c, to a tolerance relative to
# its length
function _on_segment(a, b, c; tol = 1e-6)
    ac = c - a
    L2 = dot(ac, ac)
    L2 == 0 && return false
    t = dot(b - a, ac) / L2
    return -tol <= t <= 1 + tol && norm(b - a - t * ac) <= tol * sqrt(L2)
end

# The faces of the section's cells that lie in its plane (top, bottom, lateral
# along the section) as side-on segments, each once, with their kinds
function _section_segments(mesh, ijk, bdir, columns, cells, along, sz)
    lateral = along == :i ? (Int8(1), Int8(2)) : (Int8(3), Int8(4))
    dirs = (lateral..., _TOP, _BOTTOM)
    sealed, leaky, seen = Set(columns.sealed), Set(columns.leaky), Set{Int}()
    segments, kind = NTuple{2, NTuple{2, Float64}}[], Symbol[]
    for c in cells
        for f in mesh.faces.cells_to_faces[c]
            d, o = _interior_face_direction(mesh, ijk, c, f)
            (d in dirs && !(f in seen)) || continue
            push!(seen, f)
            push!(segments, _side_on(mesh, mesh.faces.faces_to_nodes[f], d, sz))
            push!(kind, _section_face_kind(f, c, o, columns, sealed, leaky))
        end
        # boundary faces: a block's outline
        for f in mesh.boundary_faces.cells_to_faces[c]
            bdir[f] in dirs || continue
            nodes = mesh.boundary_faces.faces_to_nodes[f]
            push!(segments, _side_on(mesh, nodes, bdir[f], sz))
            push!(kind, :outline)
        end
    end
    return (segments, kind)
end

# A face side-on, from its corners: a top or bottom face split along the
# section, a lateral one by depth
_side_on(mesh, nodes, dir, sz) =
    _split_ends(sz.(_face_corners(mesh, nodes)), dir in (_TOP, _BOTTOM) ? 1 : 2)

# The kind of interior face f between cell c and cell o: a barrier, inside one
# VE block, or a block's outline
function _section_face_kind(f, c, o, columns, sealed, leaky)
    f in sealed && return :sealed
    f in leaky && return :leaky
    b = columns.partition[c]
    return b == columns.partition[o] && columns.is_ve[b] ? :internal : :outline
end

# The means of the lower and the upper half of points p, ordered by their
# coordinate `axis`
function _split_ends(p, axis)
    q = sort(p; by = x -> x[axis])
    h = length(q) ÷ 2
    mean2(v) = (sum(first, v) / length(v), sum(last, v) / length(v))
    return (mean2(q[1:h]), mean2(q[(h + 1):end]))
end

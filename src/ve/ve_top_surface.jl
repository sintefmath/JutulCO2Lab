# The top surface of a VE domain as polygons in map view, and contour lines on
# it. Geometry only; plotting is in the Makie extension.
#
# co2lab's mapPlot resamples every field onto a 200 × 200 raster and clips it
# to one boundary polygon; here fields are drawn on the columns' own top faces
# and contoured on them exactly, which also handles holes and grids in
# several pieces.

# ------------------------------- The top surface -------------------------------

"""
    ve_top_surface(mesh, columns::VEColumns; blocks = 1:number_of_columns(columns))
        -> NamedTuple

The caprock of a VE domain in map view, as one polygon per column: the top
face of its top cell. The basis for map plots and contour lines
([`ve_contour_lines`](@ref)).

# Arguments

- `mesh`, `columns`: the 3D mesh and its column partition
  ([`ve_columns`](@ref)).
- `blocks`: the columns to include; default all.

# Returns

A NamedTuple of new arrays:

- `x`, `y`, `z`: coordinates of the surface's nodes (`z`: caprock depth).
- `polygons`: `polygons[i]` holds the node indices of the top face of column
  `blocks[i]`, in the face's order (a collapsed pillar edge gives a
  triangle).
- `xc`, `yc`: the mean of each polygon's nodes.
- `outline`: node pairs `(a, b)`, `a < b`, of the edges used by only one
  polygon: the lateral boundary, holes included.
- `mesh_nodes`: for each surface node, the node of `mesh` it is.
- `blocks`: the column of each polygon, a `Vector{Int}`.

# Notes

- With all columns of a grid with several blocks per pillar (stacked
  columns, hybrid grids) the polygons of a pillar overlap; select one sheet
  with `blocks` ([`trap_surface`](@ref) does).
- Each top cell must have exactly one top face (on corner-point grids faults
  lie between pillars); an error is thrown otherwise.
"""
function ve_top_surface(mesh, columns::VEColumns;
                        blocks = 1:number_of_columns(columns))
    bdir = _boundary_directions(mesh)
    ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
    blocks = collect(Int, blocks)
    faces = [_top_face_nodes(mesh, ijk, bdir, columns.cells[columns.pos[c]])
             for c in blocks]

    # number the nodes of the surface
    index = Dict{Int, Int}()
    for f in faces, n in f
        get!(index, n, length(index) + 1)
    end
    nodes = sort(collect(keys(index)), by = n -> index[n])
    x = [mesh.node_points[n][1] for n in nodes]
    y = [mesh.node_points[n][2] for n in nodes]
    z = [mesh.node_points[n][3] for n in nodes]
    polygons = [[index[n] for n in f] for f in faces]

    xc = [sum(x[p]) / length(p) for p in polygons]
    yc = [sum(y[p]) / length(p) for p in polygons]
    return (x = x, y = y, z = z, polygons = polygons, xc = xc, yc = yc,
            outline = _outline_edges(polygons), mesh_nodes = nodes, blocks = blocks)
end

# Nodes of cell c's top face: a boundary face, or the interior face to the
# cell above (a block below another block, or below fine cells). Throws unless
# there is exactly one
_top_face_nodes(mesh, ijk, bdir, c) = _face_nodes_toward(mesh, ijk, bdir, c, _TOP)

# Nodes of cell c's face in direction `dir` (_TOP or _BOTTOM): a boundary face
# or the interior face to the cell above or below. Throws unless there is
# exactly one
function _face_nodes_toward(mesh, ijk, bdir, c, dir)
    found = [collect(mesh.boundary_faces.faces_to_nodes[f])
             for f in mesh.boundary_faces.cells_to_faces[c] if bdir[f] == dir]
    for f in mesh.faces.cells_to_faces[c]
        d, _ = _interior_face_direction(mesh, ijk, c, f)
        d == dir && push!(found, collect(mesh.faces.faces_to_nodes[f]))
    end
    name = dir == _TOP ? "top" : "bottom"
    length(found) == 1 ||
        error("cell $c ($(ijk[c])) has $(length(found)) $name faces, expected one")
    return only(found)
end

# Edges (node pairs, smaller first) that belong to exactly one polygon
function _outline_edges(polygons)
    count = Dict{Tuple{Int, Int}, Int}()
    for p in polygons, (k, a) in enumerate(p)
        b = p[k == length(p) ? 1 : k + 1]
        a == b && continue
        e = minmax(a, b)
        count[e] = get(count, e, 0) + 1
    end
    return sort!([e for (e, n) in count if n == 1])
end

# ------------------------------- Contour lines --------------------------------

"""
    ve_contour_lines(surface, values, levels) -> Vector{Vector{NTuple{4, Float64}}}

Contour lines of a field on a VE top surface, traced exactly on the columns'
own top faces.

# Arguments

- `surface`: the top surface ([`ve_top_surface`](@ref)).
- `values`: the field, one value per polygon (column) or per node of
  `surface`.
- `levels`: the contour levels.

# Returns

A new vector with, for each level in `levels`, a vector of line segments
`(x1, y1, x2, y2)`.

# Notes

- Per-column values are given to each node as the mean of the columns around
  it. Each polygon is split into triangles around its centre (the mean of its
  nodes), which takes the column's value, or for per-node values the mean of
  its nodes' values; lines are traced linearly in each triangle (marching
  triangles).
- Triangles with a `NaN` value are skipped.
"""
function ve_contour_lines(surface, values::AbstractVector, levels)
    node, centre = _node_and_centre_values(surface, values)
    x, y = surface.x, surface.y
    out = [NTuple{4, Float64}[] for _ in levels]
    # each polygon as a fan of triangles around its centre
    for (col, p) in enumerate(surface.polygons)
        xc, yc, vc = surface.xc[col], surface.yc[col], centre[col]
        for k in eachindex(p)
            a, b = p[k], p[k == length(p) ? 1 : k + 1]
            a == b && continue
            tri = ((xc, yc, vc), (x[a], y[a], node[a]), (x[b], y[b], node[b]))
            for (l, level) in enumerate(levels)
                _triangle_segment!(out[l], tri, level)
            end
        end
    end
    return out
end

# Node values and polygon-centre values of a field given per column or node
function _node_and_centre_values(surface, values)
    nn, ncol = length(surface.x), length(surface.polygons)
    if length(values) == nn
        node = collect(Float64, values)
        centre = [sum(node[p]) / length(p) for p in surface.polygons]
        return (node, centre)
    end
    length(values) == ncol ||
        throw(ArgumentError("values must be per column ($ncol) or per node ($nn)"))

    # per-column values: each node takes the mean of its columns'
    total, count = zeros(nn), zeros(Int, nn)
    for (col, p) in enumerate(surface.polygons), n in unique(p)
        total[n] += values[col]
        count[n] += 1
    end
    return (total ./ count, collect(Float64, values))
end

"Add the segment of `level` in triangle `tri`, if any. Modifies: `segments`."
function _triangle_segment!(segments, tri, level)
    any(t -> isnan(t[3]), tri) && return segments
    pts = NTuple{2, Float64}[]
    for (i, j) in ((1, 2), (2, 3), (3, 1))
        (xi, yi, vi), (xj, yj, vj) = tri[i], tri[j]
        # the edge crosses the level: one end below, the other not
        (vi < level) == (vj < level) && continue
        t = (level - vi) / (vj - vi)
        push!(pts, (xi + t * (xj - xi), yi + t * (yj - yi)))
    end
    length(pts) == 2 && push!(segments, (pts[1]..., pts[2]...))
    return segments
end

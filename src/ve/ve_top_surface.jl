# The top surface of a VE domain as polygons in map view, and contour lines on
# it (roadmap Stage 5 item 5). Geometry only — plotting is in the Makie
# extension (ext/JutulCO2LabMakieExt).
#
# co2lab's mapPlot resamples every field onto a 200 × 200 raster with
# TriScatteredInterp and clips it to one boundary polygon; here fields are
# drawn on the columns' own top faces and contoured on them exactly, which
# also handles holes and grids in several pieces.

# ------------------------------- The top surface -------------------------------

"""
    ve_top_surface(mesh, columns::VEColumns; blocks = 1:number_of_columns(columns))
        -> NamedTuple

The caprock of a VE domain in map view: the top face of each column in
`blocks` (the top face of its top cell) as a polygon. With all columns of a
grid with several blocks per pillar (stacked columns, hybrid grids) the
polygons of a pillar overlap; [`trap_surface`](@ref) takes one sheet. A top
cell must have one top face (corner-point grids: faults lie between pillars);
an error is thrown otherwise.

Returns the NamedTuple `(x, y, z, polygons, xc, yc, outline, mesh_nodes,
blocks)`, all new: `x`, `y`, `z` the coordinates of the surface's nodes (`z`
the caprock depth there); `polygons[i]` the node indices of the top face of
column `blocks[i]`, in the face's order (a collapsed pillar edge gives a
triangle); `xc`, `yc` the mean of each polygon's nodes; `outline` the node
pairs of the edges used by only one polygon (the lateral boundary, holes
included); `mesh_nodes` the node of `mesh` each surface node is; `blocks`
(a `Vector{Int}`) the column of each polygon. `mesh` and `columns` are only
read.
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
# cell above (a block below another block, or below fine cells); exactly one
function _top_face_nodes(mesh, ijk, bdir, c)
    tops = [collect(mesh.boundary_faces.faces_to_nodes[f])
            for f in mesh.boundary_faces.cells_to_faces[c] if bdir[f] == _TOP]
    for f in mesh.faces.cells_to_faces[c]
        d, _ = _interior_face_direction(mesh, ijk, c, f)
        d == _TOP && push!(tops, collect(mesh.faces.faces_to_nodes[f]))
    end
    length(tops) == 1 ||
        error("cell $c ($(ijk[c])) has $(length(tops)) top faces, expected one")
    return only(tops)
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

Contour lines of a field on a [`ve_top_surface`](@ref): `values` per column
(then each node takes the mean of the columns around it) or per node of the
surface. Each polygon is split into triangles around its node mean, which
takes the column's value (or the mean of the node values), and the lines are
traced linearly in each triangle (marching triangles). Triangles with a `NaN`
value are skipped.

Returns, for each level in `levels`, a new vector of line segments
`(x1, y1, x2, y2)`. The arguments are only read.
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

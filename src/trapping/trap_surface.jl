# The top surface as the trapping analysis sees it (roadmap Stage 8 item 1) —
# what co2lab-spillpoint takes from the top-surface grid Gt: node coordinates
# and depths, edges with the columns on either side, and per column its
# caprock depth, height, area, centroid and corner nodes.
#
# Nodes are numbered as co2lab numbers them (sorted by x, then y, then depth),
# so that ties between equally steep neighbours are broken the same way.

# ---------------------------------- Surface -----------------------------------

"""
    trap_surface(mesh, columns::VEColumns) -> NamedTuple

The caprock of a VE domain for trap analysis: [`ve_top_surface`](@ref)'s
polygons with co2lab's numbering and the quantities co2lab-spillpoint reads
from its top-surface grid `Gt`.

Returns the NamedTuple `(x, y, z, polygons, corners, edges, edge_columns, zc,
H, area, cx, cy, outline)`, all new:
- `x`, `y`, `z`: node coordinates and caprock depths, nodes sorted by
  `(x, y, z)` as co2lab's `uniqueNodes`;
- `polygons[col]`: the column's top-face nodes;
- `corners`: `4 × ncol`, each column's west-south, east-south, west-north and
  east-north nodes (co2lab's `cellNodes`; for columns with four distinct
  top-face nodes, else 0);
- `edges`: `2 × nedge` node pairs (smaller first) in co2lab's order — when
  every column has four corners, as its `cartesian_face_reordering`: each
  column's west and east edges column by column, then all south edges,
  then all north edges; otherwise sorted — and
  `edge_columns`: `2 × nedge`, the columns on either side (0 at the
  boundary);
- `zc`, `H`: caprock depth and height per column
  ([`ve_column_geometry`](@ref));
- `area`, `cx`, `cy`: the area and centroid of each column's polygon in map
  view (co2lab's `Gt.cells.volumes`, `Gt.cells.centroids`);
- `outline`: boundary edges as node pairs, as [`ve_top_surface`](@ref).
`mesh` and `columns` are only read.
"""
function trap_surface(mesh, columns::VEColumns)
    s = ve_top_surface(mesh, columns)
    geo = ve_column_geometry(mesh, columns)

    # co2lab's node order: sorted by x, y, z
    order = sortperm(collect(zip(s.x, s.y, s.z)))
    new = invperm(order)
    polygons = [new[p] for p in s.polygons]
    x, y, z = s.x[order], s.y[order], s.z[order]

    edges, edge_columns = _surface_edges(polygons)
    corners = _polygon_corners(x, y, polygons)
    # co2lab's order on lattice-like surfaces (ties in the node-based method)
    if all(>(0), corners)
        order = _co2lab_edge_order(edges, corners)
        edges, edge_columns = edges[:, order], edge_columns[:, order]
    end
    area, cx, cy = _polygon_geometry(x, y, polygons)
    outline = sort!([minmax(new[a], new[b]) for (a, b) in s.outline])
    return (x = x, y = y, z = z, polygons = polygons,
            corners = corners, edges = edges,
            edge_columns = edge_columns, zc = geo.z_top, H = geo.H,
            area = area, cx = cx, cy = cy, outline = outline)
end

# ------------------------------- Edges, corners -------------------------------

# Edges as sorted node pairs (2 × ne) and the columns on either side (2 × ne,
# the smaller column first; 0 for a boundary edge)
function _surface_edges(polygons)
    sides = Dict{Tuple{Int, Int}, Vector{Int}}()
    for (col, p) in enumerate(polygons), (k, a) in enumerate(p)
        b = p[k == length(p) ? 1 : k + 1]
        a == b && continue
        push!(get!(sides, minmax(a, b), Int[]), col)
    end
    pairs = sort!(collect(keys(sides)))
    edges = [e[i] for i in 1:2, e in pairs]
    cols = [i <= length(sides[e]) ? sort(sides[e], rev = true)[i] : 0
            for i in 1:2, e in pairs]
    # boundary edges: (column, 0); interior: (smaller, larger)
    for k in axes(cols, 2)
        cols[2, k] != 0 && (cols[:, k] = sort(cols[:, k]))
    end
    return (edges, cols)
end

# co2lab's edge order on a lattice (topSurfaceGrid's cartesian_face_reordering):
# west and east edge of each column in turn, then the south edges, then the
# north edges, each edge where it first appears
function _co2lab_edge_order(edges, corners)
    index = Dict((edges[1, k], edges[2, k]) => k for k in axes(edges, 2))
    e(a, b) = index[minmax(a, b)]
    ws, es, wn, en = (corners[i, :] for i in 1:4)
    first = [k for c in axes(corners, 2) for k in (e(ws[c], wn[c]), e(es[c], en[c]))]
    second = [[e(ws[c], es[c]) for c in axes(corners, 2)];
              [e(wn[c], en[c]) for c in axes(corners, 2)]]
    order = unique([unique(first); unique(second)])
    length(order) == size(edges, 2) || return collect(axes(edges, 2))
    return order
end

# West-south, east-south, west-north, east-north node of each column
# (4 × ncol; 0 where the polygon does not have four distinct nodes)
function _polygon_corners(x, y, polygons)
    corners = zeros(Int, 4, length(polygons))
    for (col, p) in enumerate(polygons)
        q = unique(p)
        length(q) == 4 || continue
        # the two southern nodes, then the two northern, each west to east
        s = sort(q, by = n -> (y[n], x[n]))
        south, north = sort(s[1:2], by = n -> x[n]), sort(s[3:4], by = n -> x[n])
        corners[:, col] = [south; north]
    end
    return corners
end

# Area and centroid of each polygon in map view (shoelace formula, relative to
# the polygon's first node: map coordinates of ~1e6 m would cancel badly)
function _polygon_geometry(x, y, polygons)
    n = length(polygons)
    area, cx, cy = zeros(n), zeros(n), zeros(n)
    for (col, p) in enumerate(polygons)
        x0, y0 = x[p[1]], y[p[1]]
        a, sx, sy = 0.0, 0.0, 0.0
        for (k, i) in enumerate(p)
            j = p[k == length(p) ? 1 : k + 1]
            xi, yi, xj, yj = x[i] - x0, y[i] - y0, x[j] - x0, y[j] - y0
            w = xi * yj - xj * yi
            a += w
            sx += (xi + xj) * w
            sy += (yi + yj) * w
        end
        area[col], cx[col], cy[col] = abs(a) / 2, x0 + sx / (3a), y0 + sy / (3a)
    end
    return (area, cx, cy)
end

# Sheets of a VE grid: the caprock surfaces of grids with several VE blocks
# per pillar (stacked columns below holes, sealed or leaky layers, VE regions
# under fine cells). A sheet is a contiguous surface of connected top faces of
# VE columns; each gets its own trap analysis (trap_surface(...; sheet)), and
# CO2 leaving a sheet over an open rim edge spills up to the sheet above
# (sheet_spill_targets).

# ---------------------------------- Sheets ------------------------------------

"""
    ve_sheets(mesh, columns::VEColumns) -> NamedTuple

The sheets of a VE grid: two VE columns in different pillars belong to the
same sheet when their top cells share a face of `mesh` (their top faces join,
or meet across a fault with partial offset); a sheet is a connected set of
such columns. Columns lower in a pillar (below a hole, a sealed or leaky
face, or fine cells) form sheets of their own; fine cells (`is_ve` false)
belong to none. Sheets are numbered by their first column (columns are in
pillar order, top first).

Returns the NamedTuple `(sheet, nsheets, top)`, all new: `sheet` the sheet of
each column (`Vector{Int}`, 0 for fine cells), the number of sheets, and per
sheet whether all its columns are the first block of their pillar
(`Vector{Bool}`: the formation's top, nothing above it). `mesh` and
`columns` are only read.
"""
function ve_sheets(mesh, columns::VEColumns)
    nb = number_of_columns(columns)
    # the VE column whose top cell each cell is (0: none)
    top_of = zeros(Int, number_of_cells(mesh))
    for b in 1:nb
        columns.is_ve[b] && (top_of[columns.cells[columns.pos[b]]] = b)
    end

    # top cells sharing a face, in different pillars: one sheet
    g = Graphs.SimpleGraph(nb)
    for (c1, c2) in eachcol(get_neighborship(mesh))
        a, b = top_of[c1], top_of[c2]
        a > 0 && b > 0 && columns.column_ij[a] != columns.column_ij[b] &&
            Graphs.add_edge!(g, a, b)
    end
    comps = sort!([sort(c) for c in Graphs.connected_components(g)
                   if columns.is_ve[first(c)]], by = first)

    sheet = zeros(Int, nb)
    foreach(((k, c),) -> sheet[c] .= k, enumerate(comps))
    first_block = _first_blocks(columns)
    top = [all(b -> first_block[b], c) for c in comps]
    return (sheet = sheet, nsheets = length(comps), top = top)
end

# Per column: whether it is the first block of its pillar and VE (nothing
# above it); columns of a pillar are consecutive, top first
function _first_blocks(columns::VEColumns)
    ij = columns.column_ij
    return [columns.is_ve[b] && (b == 1 || ij[b-1] != ij[b])
            for b in 1:number_of_columns(columns)]
end

# The columns of the trap surface `sheet`: `:top` the first block of each
# pillar if VE (the formation's top: on grids with one block per pillar, all
# columns), an `Integer` the columns of that sheet of `ve_sheets`
function _surface_blocks(mesh, columns::VEColumns, sheet)
    sheet === :top && return findall(_first_blocks(columns))
    sheet isa Integer ||
        throw(ArgumentError("sheet must be :top or a sheet number, got $sheet"))
    sh = ve_sheets(mesh, columns)
    1 <= sheet <= sh.nsheets ||
        throw(ArgumentError("sheet $sheet: the grid has $(sh.nsheets) sheets"))
    return findall(==(sheet), sh.sheet)
end

# ------------------------------ Spill between sheets ---------------------------

"""
    sheet_spill_targets(mesh, columns::VEColumns, s) -> Vector{NamedTuple}

Where CO2 goes that reaches an open rim edge of the trap surface `s`
([`trap_surface`](@ref) of `mesh` and `columns`, typically one sheet): across
the edge into the cell of the neighbouring pillar beside the top of the
polygon's column (the shallowest one touching it), then up that pillar to the
first VE column — the CO2 rises in it to its top, on the sheet that column
belongs to. Sealed edges (faults) are closed and not listed.

Returns a new vector with one NamedTuple `(edge, polygon, block, target,
target_sheet)` per open boundary edge of `s`: the edge (index into
`s.edges`), the polygon of `s` it bounds and that polygon's column, the
column the CO2 reaches (`target`; 0 if the edge is on the grid's outer
boundary, where CO2 leaves the domain, or the pillar beyond has only fine
cells up to its top, where it stays in the fine region: then the topmost fine
cell's column) and its sheet ([`ve_sheets`](@ref); 0 for none). The edge's
side is the lateral side of the top cell whose faces hold both its end
nodes. `mesh`, `columns` and `s` are only read.
"""
function sheet_spill_targets(mesh, columns::VEColumns, s)
    sh = ve_sheets(mesh, columns)
    ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
    bdir = _boundary_directions(mesh)
    depth = tpfv_geometry(mesh).cell_centroids[3, :]
    T = NamedTuple{(:edge, :polygon, :block, :target, :target_sheet),
                   NTuple{5, Int}}
    out = T[]
    for k in axes(s.edges, 2)
        # boundary edges only, and only open ones
        s.edge_columns[2, k] == 0 && !_is_sealed(s, k) || continue
        p = s.edge_columns[1, k]
        b = s.blocks[p]
        top = columns.cells[columns.pos[b]]
        ends = [(s.x[n], s.y[n], s.z[n]) for n in s.edges[:, k]]
        d = _edge_direction(mesh, ijk, bdir, top, ends)
        target = d == 0 ? 0 : _spill_target(mesh, columns, ijk, depth, b, d)
        tsheet = target == 0 ? 0 : sh.sheet[target]
        push!(out, T((k, p, b, target, tsheet)))
    end
    return out
end

# The lateral side (1..4: left, right, lower, upper) of cell `top` whose
# faces hold both points `ends` (x, y, z) as nodes; 0 if none
function _edge_direction(mesh, ijk, bdir, top, ends)
    pos(n) = Tuple(mesh.node_points[n])
    sides = _lateral_face_nodes(mesh, ijk, bdir, top)
    for d in 1:4
        pts = Set(pos(n) for (dd, nodes) in sides if dd == d for n in nodes)
        all(in(pts), ends) && return d
    end
    return 0
end

# The column CO2 reaches from the top of column b over its side d: the
# shallowest cell beyond that side's faces, then up its pillar to the first
# VE column (0: no cell beyond; the topmost fine cell's column if all fine)
function _spill_target(mesh, columns, ijk, depth, b, d)
    top = columns.cells[columns.pos[b]]
    beyond = Int[]
    for f in mesh.faces.cells_to_faces[top]
        dir, o = _interior_face_direction(mesh, ijk, top, f)
        dir == d && push!(beyond, o)
    end
    isempty(beyond) && return 0
    t = columns.partition[beyond[argmin(depth[beyond])]]
    t == 0 && return 0
    # up the pillar: the columns above are the preceding ones of that pillar
    while !columns.is_ve[t] && t > 1 && columns.column_ij[t-1] == columns.column_ij[t]
        t -= 1
    end
    return t
end

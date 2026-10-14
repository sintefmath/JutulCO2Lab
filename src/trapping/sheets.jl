# Sheets: the separate caprock surfaces of a VE grid with several VE blocks
# per pillar (stacked columns below holes, sealed or leaky layers, VE regions
# under fine cells). Each sheet gets its own trap analysis
# (`trap_surface(...; sheet)`); CO2 spilling over an open rim edge of a sheet
# rises to the sheet above (`sheet_spill_targets`).

# ---------------------------------- Sheets ------------------------------------

"""
    ve_sheets(mesh, columns::VEColumns) -> NamedTuple

Split the VE columns into sheets: connected caprock surfaces, each analysed
for traps on its own (`trap_surface(...; sheet)`).

Two VE columns in different pillars are on the same sheet when their top
cells share a face of `mesh` (top faces that join, or meet across a fault
with partial offset). Columns lower in a pillar (below a hole, a sealed or
leaky face, or fine cells) form sheets of their own; fine cells belong to
none.

# Arguments

- `mesh`: the 3D mesh.
- `columns`: its column partition ([`ve_columns`](@ref)).

# Returns

A new NamedTuple:
- `sheet`: the sheet of each column (`Vector{Int}`; 0 for fine cells);
- `nsheets`: the number of sheets;
- `top`: per sheet, whether all its columns are the first block of their
  pillar, i.e. the formation's top with nothing above (`Vector{Bool}`).

Sheets are numbered by their first column (columns are in pillar order, top
first).
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

# The columns of the trap surface `sheet`: `:top` gives the first block of each
# pillar if VE (the formation's top; all columns on grids with one block per
# pillar), an `Integer` the columns of that sheet of `ve_sheets`
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

Where CO2 goes when it reaches an open rim edge of the trap surface `s`
(typically one sheet): across the edge into the neighbouring pillar, then up
that pillar to the first VE column, in which it rises to the top of that
column's sheet. This links the trap analyses of stacked sheets.

# Arguments

- `mesh`, `columns`: the 3D mesh and its column partition.
- `s`: a [`trap_surface`](@ref) of `mesh` and `columns`.

# Returns

A new vector with one NamedTuple per open boundary edge of `s` (sealed
edges, e.g. faults, are closed and not listed):
- `edge`: the edge (index into `s.edges`);
- `polygon`: the polygon of `s` it bounds;
- `block`: that polygon's column;
- `target`: the column the CO2 reaches; 0 if the edge is on the grid's outer
  boundary (CO2 leaves the domain). If the pillar beyond has only fine cells
  up to its top, the CO2 stays in the fine region: the topmost fine cell's
  column;
- `target_sheet`: the sheet of `target` ([`ve_sheets`](@ref)); 0 for none.

# Notes

- The CO2 enters the neighbouring pillar at its shallowest cell beside the
  top cell of the edge's column. The edge's side is the lateral side of that
  top cell whose faces hold both end nodes of the edge.
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

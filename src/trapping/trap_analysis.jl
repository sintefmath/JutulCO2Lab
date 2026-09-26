# The trapping analysis's front end (roadmap Stage 8 item 4): co2lab's
# trapAnalysis(Gt, method, ...), on a trap surface.

"""
    trap_analysis(s; method = :cell, closed_boundary_edges = [],
                  closed_fault_edges = [], project_to_cells = true)
        -> NamedTuple

Structural traps of a [`trap_surface`](@ref) `s`, as co2lab's
`trapAnalysis(Gt, method == :cell, ...)`:

- `method = :cell`: from column depths and column adjacency (co2lab's
  cell-based method); `method = :node`: along the surface's nodes and edges
  (the node-based method, co2lab's default in its VE examples).
- `closed_boundary_edges`, `closed_fault_edges`: edges (node pairs of
  `s.edges`) across which CO2 does not spill — node-based method only; the
  cell-based method throws if given any (co2lab only warns and ignores
  them).
- `project_to_cells`: node-based method only (`false`: results per node —
  trap, spill region per node, rivers as node lists, `top` the top nodes).

Returns a new NamedTuple `(traps, trap_z, trap_regions, trap_adj, rivers,
top)`: the trap index of each column (0: none), the spill depth of each trap,
the trap each column drains into (0: none), which trap spills into which (a
sparse `Bool` matrix, `trap_adj[i, j]`: `i` spills into `j`), the rivers
leaving each trap (`rivers[t]`: vectors of columns, from the trap
downstream) and the columns at the top of the traps. Trap numbers may differ
from co2lab's (it numbers equal-depth pieces in `dmperm`'s order). `s` is
only read.
"""
function trap_analysis(s; method::Symbol = :cell, closed_boundary_edges = [],
                       closed_fault_edges = [], project_to_cells::Bool = true)
    if method == :cell
        isempty(closed_boundary_edges) && isempty(closed_fault_edges) ||
            throw(ArgumentError("the cell-based method does not support closed " *
                                "boundary or fault edges: use method = :node"))
        return _cell_trap_analysis(s)
    end
    method == :node ||
        throw(ArgumentError("method must be :cell or :node, got :$method"))
    nt = _node_traps(s, _closed_nodes(closed_boundary_edges),
                     _closed_nodes(closed_fault_edges))
    project_to_cells && return _node_to_cell_traps(s, nt)
    top = findall(n -> nt.down[n] == 0 && nt.trap_regions[n] != 0, eachindex(s.z))
    return (traps = nt.traps, trap_z = nt.trap_z, trap_regions = nt.trap_regions,
            trap_adj = nt.trap_adj, rivers = nt.rivers, top = top)
end

# Nodes closed by a set of closed edges (node pairs): those in more than one
# of them (co2lab: the duplicates), so that the chain's end nodes stay open
function _closed_nodes(edges)
    nodes = sort!([n for e in edges for n in e])
    return unique(nodes[findall(iszero, diff(nodes))])
end

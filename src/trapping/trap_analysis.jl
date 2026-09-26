# The trapping analysis's front end (roadmap Stage 9 item 4): co2lab's
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
- `project_to_cells`: node-based method only (`false`: results per node).

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
    method == :node && throw(ArgumentError("method = :node is not implemented yet"))
    throw(ArgumentError("method must be :cell or :node, got :$method"))
end

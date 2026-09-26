# Comparisons of the trapping analysis with co2lab-spillpoint (roadmap Stage 9),
# used by the package tests (runtests.jl) and the Stage 9 gate. Each takes the
# MRST export traps.mat (scripts/mrst_dump_traps.m, read with matread) and
# returns an OracleReport.

# The mesh and columns of one exported case: an Atlas formation as co2lab builds
# it (co2lab_compatible, no thin-cell removal), an IGEMS surface as co2lab
# builds it ("IGEMS/<family>/<i>", transposed), or the synthetic sine surface
function _trap_case_mesh(e)
    name, c = e["name"], Int(e["coarsening"])
    name == "sine" && return _sine_mesh()
    if startswith(name, "IGEMS/")
        _, family, i = split(name, "/")
        m = igems_mesh(family, parse(Int, i); coarsening = c, co2lab_compatible = true)
        return (m.mesh, m.columns)
    end
    m = atlas_mesh(name; coarsening = c, co2lab_compatible = true, min_thickness = 0.0)
    return (m.mesh, m.columns)
end

# mrst_dump_traps.m's synthetic case: 30 x 20 cells of 100 m, 50 m thick, the
# caprock at 1000 - 20 sin(x/400) cos(y/300) + 0.005 x (cartGrid, nodes moved)
function _sine_mesh()
    nx, ny = 30, 20
    f(x, y) = 1000 - 20 * sin(x / 400) * cos(y / 300) + 0.005 * x
    coord = Float64[]
    for j in 0:ny, i in 0:nx
        x, y = 100.0 * i, 100.0 * j
        append!(coord, [x, y, f(x, y), x, y, f(x, y) + 50])
    end
    # each cell's eight corners: corner a of 2n along x lies on pillar 1 + a ÷ 2
    zcorn = [f(100.0 * (a ÷ 2), 100.0 * (b ÷ 2)) + 50.0 * (c - 1)
             for c in 1:2 for b in 1:2ny for a in 1:2nx]
    grid = Dict("cartDims" => (nx, ny, 1), "COORD" => coord, "ZCORN" => zcorn,
                "ACTNUM" => fill(true, nx, ny, 1))
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    return (mesh, first(ve_columns(mesh)))
end

# Item 1: the trap surface against co2lab's Gt — nodes (order included),
# edges and their columns, caprock depths, heights, areas, centroids, corners
function _trap_surface_vs_mrst(d)
    rep = OracleReport("Trap surface vs co2lab's Gt"; mrst = d)
    exact = Tolerance(0.0, 0.0)
    for e in d["cases"]
        label = "$(e["name"]) c$(Int(e["coarsening"]))"
        mesh, cols = _trap_case_mesh(e)
        s = trap_surface(mesh, cols)

        nodes = hcat(s.x, s.y, s.z)
        check!(rep, "$label: nodes (co2lab's order)", nodes,
               hcat(e["node_coords"], vec(e["node_z"])); tol = Tolerance(0.0, 1e-9))
        E = sort(Int.(e["edge_nodes"]), dims = 2)
        N = Int.(e["edge_neighbors"])
        mrst_edges = Dict(Tuple(E[k, :]) => sort(N[k, :]) for k in axes(E, 1))
        ours = Dict(Tuple(s.edges[:, k]) => sort(s.edge_columns[:, k])
                    for k in axes(s.edges, 2))
        check!(rep, "$label: edges and their columns", ours == mrst_edges ? 0 : 1, 0;
               tol = exact)
        same_corners = s.corners == permutedims(Int.(e["cell_corners"]))
        check!(rep, "$label: corner nodes", same_corners ? 0 : 1, 0; tol = exact)

        rel = Tolerance(1e-12, 1e-9)
        check!(rep, "$label: caprock depth", s.zc, vec(e["cell_z"]); tol = rel)
        check!(rep, "$label: height", s.H, vec(e["cell_H"]); tol = rel)
        check!(rep, "$label: area", s.area, vec(e["cell_areas"]); tol = rel)
        check!(rep, "$label: centroid", hcat(s.cx, s.cy), e["cell_centroids"]; tol = rel)
    end
    return rep
end

# The trap surface of one exported case with co2lab's column depths: our depths
# agree with co2lab's only to about an ulp (face centroids computed in a
# different order), and on surfaces with rounded depths (IGEMS: 1 cm) exact
# ties are common, so an ulp can decide a trap. The algorithms are therefore
# compared on identical input; _trap_surface_vs_mrst checks the depths.
function _trap_case_surface(e)
    mesh, cols = _trap_case_mesh(e)
    return (; trap_surface(mesh, cols)..., zc = _vecs(e["cell_z"]))
end

# Our trap numbers → co2lab's, matching traps by their columns (co2lab's trap
# numbering follows dmperm and is arbitrary); nothing if the sets differ
function _trap_mapping(ours, theirs)
    n = maximum(ours; init = 0)
    m = zeros(Int, n)
    for t in 1:n
        ids = unique(theirs[ours .== t])
        length(ids) == 1 && ids[1] > 0 || return nothing
        m[t] = ids[1]
    end
    allunique(m) && count(>(0), theirs) == count(>(0), ours) || return nothing
    return m
end

# co2lab's rivers per trap (exported as [trap; river; cells...]) as sets
_mrst_rivers(t) = [(Int(r[1]), Int.(r[3:end])) for r in vec(t["rivers"]) if !isempty(r)]

# Item 2 (cell-based) or 3 (node-based): trapAnalysis's result on identical
# input (_trap_case_surface), traps matched by their columns — trap
# membership, spill depths, spill regions, trap connections, tops and rivers
function _traps_vs_mrst(d, method::Symbol)
    rep = OracleReport("Traps ($method-based) vs co2lab's trapAnalysis"; mrst = d)
    exact = Tolerance(0.0, 0.0)
    key = method == :cell ? "cell_based" : "node_based"
    for e in d["cases"]
        label = "$(e["name"]) c$(Int(e["coarsening"]))"
        r = trap_analysis(_trap_case_surface(e); method)
        t = e[key]
        theirs = Int.(vec(t["traps"]))
        check!(rep, "$label: number of traps", maximum(r.traps; init = 0),
               maximum(theirs; init = 0); tol = exact)
        m = _trap_mapping(r.traps, theirs)
        check!(rep, "$label: trap columns", isnothing(m) ? 1 : 0, 0; tol = exact)
        isnothing(m) && continue

        check!(rep, "$label: spill depths", r.trap_z, vec(t["trap_z"])[m];
               tol = Tolerance(1e-12, 0.0))
        regions = [g == 0 ? 0 : m[g] for g in r.trap_regions]
        check!(rep, "$label: spill regions",
               count(regions .!= Int.(vec(t["trap_regions"]))), 0; tol = exact)
        adj = Set((m[i], m[j]) for (i, j, v) in zip(findnz(r.trap_adj)...) if v)
        A = Int.(t["trap_adj"])
        same_adj = adj == Set(Tuple(A[k, :]) for k in axes(A, 1))
        check!(rep, "$label: trap connections", same_adj ? 0 : 1, 0; tol = exact)
        check!(rep, "$label: tops", sort(r.top) == sort(Int.(vec(t["top"]))) ? 0 : 1, 0;
               tol = exact)
        ours = Set((m[i], cells) for i in eachindex(r.rivers) for cells in r.rivers[i])
        check!(rep, "$label: rivers", ours == Set(_mrst_rivers(t)) ? 0 : 1, 0;
               tol = exact)
    end
    return rep
end

# A MATLAB vector as read by matread: scalars for one element, matrices for more
_vecs(x) = x isa Number ? [x] : vec(x)

# co2lab's trees as (root, traps set, value), mapped through m⁻¹ = identity on
# co2lab's numbers
_mrst_trees(T) = [(Int(r), Set(Int.(_vecs(t))), v) for (r, t, v) in
                  zip(_vecs(T["root"]), _vecs(T["traps"]), _vecs(T["value"]))]

# Item 5: volumes, downstream traps, capacity trees, injection point and the
# flattened surface, on the traps of `method`, mapped to co2lab's numbers
function _trap_quantities_vs_mrst(d, method::Symbol)
    rep = OracleReport("Trap quantities ($method-based) vs co2lab"; mrst = d)
    exact, rel = Tolerance(0.0, 0.0), Tolerance(1e-12, 0.0)
    # sums over up to thousands of columns (measured 3e-15)
    vol = Tolerance(1e-12, 0.0)
    for e in d["cases"]
        label = "$(e["name"]) c$(Int(e["coarsening"]))"
        t = e[method == :cell ? "cell_based" : "node_based"]
        isempty(t["volumes"]) && continue
        s = _trap_case_surface(e)
        r = trap_analysis(s; method)
        m = _trap_mapping(r.traps, Int.(_vecs(t["traps"])))
        if isnothing(m)
            check!(rep, "$label: trap columns", 1, 0; tol = exact)
            continue
        end

        check!(rep, "$label: volumes", trap_volumes(s, r), _vecs(t["volumes"])[m];
               tol = vol)
        down = [Set(m[downstream_traps(r, i)]) for i in eachindex(m)]
        theirs = [Set(Int.(_vecs(t["downstream"][m[i]]))) for i in eachindex(m)]
        check!(rep, "$label: downstream traps", count(down .!= theirs), 0; tol = exact)
        for (key, kw) in (("trees", (;)), ("trees_all", (remove_overlap = false,
                                                         calculate_all = true)))
            ours = maximize_trapping(s, r; kw...)
            mt = _mrst_trees(t[key])
            check!(rep, "$label: $key roots and traps",
                   [(m[x.root], Set(m[x.traps])) for x in ours] ==
                   [(x[1], x[2]) for x in mt] ? 0 : 1, 0; tol = exact)
            check!(rep, "$label: $key values", [x.value for x in ours], last.(mt);
                   tol = vol)
        end
        inj = optimal_injection_point(s, r)
        check!(rep, "$label: injection column", inj.column, t["injection"]["cell"];
               tol = exact)
        check!(rep, "$label: injection volume", inj.volume, t["injection"]["volume"];
               tol = vol)
        f = flatten_traps(s, r)
        check!(rep, "$label: flattened surface", [f.zc; f.z],
               [_vecs(t["flat_cell_z"]); _vecs(t["flat_node_z"])]; tol = rel)
    end
    return rep
end

# Closed boundary and fault edges (node-based): the sine surface with its west
# boundary and a fault line closed, against co2lab (node numbers agree)
function _closed_edges_vs_mrst(d)
    rep = OracleReport("Closed edges (node-based) vs co2lab"; mrst = d)
    exact = Tolerance(0.0, 0.0)
    c = d["closed"]
    s = trap_surface(_sine_mesh()...)
    pairs(M) = [Tuple(sort(Int.(M[k, :]))) for k in axes(M, 1)]
    r = trap_analysis(s; method = :node,
                      closed_boundary_edges = pairs(c["closed_boundary"]),
                      closed_fault_edges = pairs(c["closed_fault"]))
    open = trap_analysis(s; method = :node)
    m = _trap_mapping(r.traps, Int.(_vecs(c["traps"])))
    check!(rep, "trap columns", isnothing(m) ? 1 : 0, 0; tol = exact)
    isnothing(m) && return rep
    check!(rep, "spill depths", r.trap_z, _vecs(c["trap_z"])[m];
           tol = Tolerance(1e-12, 0.0))
    regions = [g == 0 ? 0 : m[g] for g in r.trap_regions]
    check!(rep, "spill regions", count(regions .!= Int.(_vecs(c["trap_regions"]))), 0;
           tol = exact)
    # the closed edges change the result (the check would be empty otherwise)
    same = r.trap_regions == open.trap_regions
    check!(rep, "differs from the open case", same ? 1 : 0, 0; tol = exact)
    return rep
end

# Comparisons of the trapping analysis with co2lab-spillpoint, used by the
# package tests (runtests.jl). Each takes the MRST export traps.mat and returns
# an OracleReport.

# The mesh and columns of one exported case: an Atlas formation as co2lab builds
# it (co2lab_compatible, no thin-cell removal), an IGEMS surface as co2lab
# builds it ("IGEMS/<family>/<i>", transposed), a synthetic surface, or one of
# the grids of co2lab-spillpoint's examples (the export's 'examples')
function _trap_case_mesh(e)
    name, c = e["name"], Int(e["coarsening"])
    name == "sine" && return _sine_mesh()
    name == "shoebox" && return _shoebox_mesh()
    name == "synthetic101" && return _synthetic101_mesh()
    if startswith(name, "IGEMS/")
        _, family, i = split(name, "/")
        m = igems_mesh(family, parse(Int, i); coarsening = c, co2lab_compatible = true)
        return (m.mesh, m.columns)
    elseif startswith(name, "Johansen/")
        m = johansen_mesh(Symbol(split(name, "/")[2]))
        return (m.mesh, m.columns)
    elseif startswith(name, "sleipner")
        grid = sleipner_grid_section(; min_perm = 0.0)
        mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
        return (mesh, first(ve_columns(mesh)))
    elseif name == "Huginfmwest/south"
        return _hugin_west_south(c)
    elseif name == "Utsirafm/thick"
        m = atlas_mesh("Utsirafm"; coarsening = c, co2lab_compatible = true)
        return (m.mesh, m.columns)
    end
    m = atlas_mesh(name; coarsening = c, co2lab_compatible = true, min_thickness = 0.0)
    return (m.mesh, m.columns)
end

# A faulted synthetic grid: 10 x 7 cells of 100 m, 20 m thick, an anticline
# rising east (crest along y = 350 m); a fault between cells 5 and 6 with the
# east block `throw` metres down (40: complete offset)
function _faulted_mesh(throw)
    nx, ny, h = 10, 7, 100.0
    f(x, y) = 1000 - 0.02x + 2e-4 * (y - 350)^2
    coord = Float64[]
    for j in 0:ny, i in 0:nx
        append!(coord, [h * i, h * j, 900.0, h * i, h * j, 1200.0])
    end
    # corner a of 2nx along x: pillar 1 + a ÷ 2, cell 1 + (a - 1) ÷ 2
    top(a, b) = f(h * (a ÷ 2), h * (b ÷ 2)) + ((a + 1) ÷ 2 > nx ÷ 2 ? throw : 0.0)
    zcorn = [top(a, b) + 20.0 * (c - 1) for c in 1:2 for b in 1:2ny for a in 1:2nx]
    grid = Dict("cartDims" => (nx, ny, 1), "COORD" => coord, "ZCORN" => zcorn,
                "ACTNUM" => fill(true, nx, ny, 1))
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    return (mesh, first(ve_columns(mesh)))
end

# trapsHuginWest's subset: cells with centroids east of 435 km, south of 6450 km
function _hugin_west_south(c)
    m = atlas_mesh("Huginfmwest"; coarsening = c, co2lab_compatible = true,
                   min_thickness = 0.0)
    x = tpfv_geometry(m.mesh).cell_centroids
    sub = grid_section_subset(m.grid, m.mesh, (x[1, :] .> 4.35e5) .& (x[2, :] .< 6.45e6))
    return (sub.mesh, sub.columns)
end

# firstTrappingExample's sandbox: 100 x 100 cells over 10 km x 5 km, 50 m
# thick at 1000 m, dipping in x, the top perturbed (cartGrid, nodes moved)
function _shoebox_mesh()
    x = JutulCO2Lab._linspace(0.0, 10000.0, 101)
    y = JutulCO2Lab._linspace(0.0, 5000.0, 101)
    xn, yn = x ./ 10000, y ./ 5000
    zt = [a - 0.2 * sin(5 * pi * a) * sin(5 * pi * b^1.5) - 0.075 * sin(1.25 * pi * b) +
          0.15 * sin(a + b) for a in xn, b in yn]
    zb = [1 + a for a in xn, b in yn]
    m = surface_mesh(x, y, zt .* 50 .+ 1000; bottom = zb .* 50 .+ 1000)
    return (m.mesh, m.columns)
end

# showTrapsInteractively's synthetic grid: 101 x 101 cells over 1 km, 10 m
# thick at about 100 m, the top a sloping egg-box
function _synthetic101_mesh()
    x = JutulCO2Lab._linspace(0.0, 1000.0, 102)
    f(a, b, z) = 100 + z + a * 0.01 - 2 * sin(pi * a / 200) * sin(pi * b / 200)
    m = surface_mesh(x, x, [f(a, b, 0.0) for a in x, b in x];
                     bottom = [f(a, b, 10.0) for a in x, b in x])
    return (m.mesh, m.columns)
end

# The export's synthetic case: 30 x 20 cells of 100 m, 50 m thick, the
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
    for e in d["cases"]
        label = "$(e["name"]) c$(Int(e["coarsening"]))"
        mesh, cols = _trap_case_mesh(e)
        # co2lab stitches faults (open)
        s = trap_surface(mesh, cols; faults = :open)
        ncol = length(_vecs(e["cell_z"]))
        check!(rep, "$label: number of columns", length(s.zc), ncol;
               tol = Tolerance(0.0, 0.0))
        length(s.zc) == ncol || continue
        _surface_topology_vs_mrst!(rep, label, s, e)
        _column_quantities_vs_mrst!(rep, label, s, e, _split_columns(mesh, cols))
    end
    return rep
end

# Nodes, edges with their columns, and corners of one case against co2lab's.
# Modifies: `rep` (checks added).
function _surface_topology_vs_mrst!(rep, label, s, e)
    exact = Tolerance(0.0, 0.0)
    # the detrended Sleipner surface: MRST's own depths (by design)
    own = e["name"] == "sleipner/detrended"
    nodes = own ? hcat(s.x, s.y) : hcat(s.x, s.y, s.z)
    theirs = own ? e["node_coords"] : hcat(e["node_coords"], vec(e["node_z"]))
    check!(rep, "$label: nodes (co2lab's order)", nodes, theirs;
           tol = Tolerance(0.0, 1e-9))
    E = sort(Int.(e["edge_nodes"]), dims = 2)
    N = Int.(e["edge_neighbors"])
    mrst_edges = Dict(Tuple(E[k, :]) => sort(N[k, :]) for k in axes(E, 1))
    ours = Dict(Tuple(s.edges[:, k]) => sort(s.edge_columns[:, k])
                for k in axes(s.edges, 2))
    check!(rep, "$label: edges and their columns", ours == mrst_edges ? 0 : 1, 0;
           tol = exact)
    # co2lab exports corners only if every cell has four; their order (logical
    # on co2lab's side, geometric here) differs on rotated grids
    if haskey(e, "cell_corners")
        C = permutedims(Int.(e["cell_corners"]))
        same = all(sort(s.corners[:, c]) == sort(C[:, c]) for c in axes(C, 2))
        check!(rep, "$label: corner nodes (as sets)", same ? 0 : 1, 0; tol = exact)
    else
        check!(rep, "$label: not all columns quadrilateral",
               all(>(0), s.corners) ? 1 : 0, 0; tol = exact)
    end
    return rep
end

# Depth, height, area and centroid per column against co2lab's, leaving out
# columns with fault-split faces (`split`: Jutul's centroids of the twisted
# faces, with the extra nodes, differ from MRST's by design) and folded top
# faces (crossing themselves in map view: co2lab adds the triangle areas
# unsigned, the shoelace formula nets them). Modifies: `rep` (checks added).
function _column_quantities_vs_mrst!(rep, label, s, e, split)
    rel = Tolerance(1e-12, 1e-9)
    kept, flat = .!split, .!_folded_columns(s)
    e["name"] == "sleipner/detrended" ||
        check!(rep, "$label: caprock depth ($(count(split)) fault-split columns " *
               "left out)", s.zc[kept], vec(e["cell_z"])[kept]; tol = rel)
    check!(rep, "$label: height", s.H[kept], vec(e["cell_H"])[kept]; tol = rel)
    check!(rep, "$label: area ($(count(.!flat)) folded columns left out)",
           s.area[flat], vec(e["cell_areas"])[flat]; tol = rel)
    check!(rep, "$label: centroid", hcat(s.cx, s.cy)[flat, :],
           e["cell_centroids"][flat, :]; tol = rel)
    return rep
end

# Columns with a cell that has a face of more than four nodes: faces split
# where a fault offsets the neighbouring cells (GeoEnergyIO)
function _split_columns(mesh, cols)
    big(nodes) = length(unique(nodes)) > 4
    split = fill(false, length(cols.pos) - 1)
    for col in eachindex(split), c in column_cells(cols, col)
        split[col] |= any(f -> big(mesh.faces.faces_to_nodes[f]),
                          mesh.faces.cells_to_faces[c]) ||
                      any(f -> big(mesh.boundary_faces.faces_to_nodes[f]),
                          mesh.boundary_faces.cells_to_faces[c])
    end
    return split
end

# Columns whose polygon folds over itself: the triangles from its node mean
# to its edges do not all turn the same way
function _folded_columns(s)
    map(s.polygons) do p
        q = unique(p)
        x0, y0 = sum(s.x[q]) / length(q), sum(s.y[q]) / length(q)
        turn(a, b) = (s.x[a] - x0) * (s.y[b] - y0) - (s.x[b] - x0) * (s.y[a] - y0)
        signs = [sign(turn(q[k], q[mod1(k + 1, length(q))])) for k in eachindex(q)]
        any(>(0), signs) && any(<(0), signs)
    end
end

# The trap surface of one exported case with co2lab's column depths, heights,
# areas and centroids: our depths agree with co2lab's only to about an ulp
# (face centroids computed in a different order), and on surfaces with rounded
# depths (IGEMS: 1 cm) exact ties are common, so an ulp can decide a trap; on
# faulted grids split and folded faces differ by design. The algorithms are
# therefore compared on identical input; _trap_surface_vs_mrst checks the
# surface itself.
function _trap_case_surface(e)
    s = trap_surface(_trap_case_mesh(e)...; faults = :open)
    c = e["cell_centroids"]
    s = (; s..., zc = _vecs(e["cell_z"]), H = _vecs(e["cell_H"]),
         area = _vecs(e["cell_areas"]), cx = c[:, 1], cy = c[:, 2])
    # the detrended Sleipner surface: MRST's node depths too
    e["name"] == "sleipner/detrended" && (s = (; s..., z = _vecs(e["node_z"])))
    return s
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
        t = e[key]
        # co2lab fails (node-based method, columns not quadrilateral): so must we
        if haskey(t, "error")
            fails = try
                trap_analysis(_trap_case_surface(e); method); false
            catch err
                err isa ArgumentError
            end
            check!(rep, "$label: fails as co2lab", fails ? 0 : 1, 0; tol = exact)
            continue
        end
        r = trap_analysis(_trap_case_surface(e); method)
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
        (haskey(t, "error") || isempty(t["volumes"])) && continue
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
        for (key, kw) in (("trees", (;)), ("trees_overlap", (remove_overlap = false,)),
                          ("trees_all", (remove_overlap = false, calculate_all = true)))
            haskey(t, key) || continue
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

# trapsJohansen's injection point: the column nearest to the benchmark well,
# on each Johansen grid exported with it
function _well_columns_vs_mrst(d)
    rep = OracleReport("Johansen injection columns vs co2lab"; mrst = d)
    w = johansen_well_point()
    for e in d["cases"]
        haskey(e, "well_column") || continue
        s = _trap_case_surface(e)
        col = argmin((s.cx .- w[1]) .^ 2 .+ (s.cy .- w[2]) .^ 2)
        check!(rep, "$(e["name"]) c$(Int(e["coarsening"]))", col, e["well_column"];
               tol = Tolerance(0.0, 0.0))
    end
    return rep
end

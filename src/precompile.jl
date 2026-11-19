# Precompiled workload: the trap analysis's steps on a small synthetic caprock
# (a sine surface, 12 x 8 columns), so that their first use in a session (e.g.
# the trap viewer's first formation) does not wait for compilation. No data
# files are read.

@setup_workload begin
    nx, ny = 12, 8
    f(x, y) = 1000 - 20 * sin(x / 400) * cos(y / 300) + 0.005 * x
    coord = Float64[]
    for j in 0:ny, i in 0:nx
        x, y = 100.0 * i, 100.0 * j
        append!(coord, [x, y, f(x, y), x, y, f(x, y) + 50])
    end
    # each cell's eight corners: corner a of 2n along x lies on pillar 1 + a ÷ 2
    zcorn = [f(100.0 * (a ÷ 2), 100.0 * (b ÷ 2)) + 50.0 * (c - 1)
             for c in 1:2 for b in 1:2ny for a in 1:2nx]
    grid = Dict{String, Any}("cartDims" => (nx, ny, 1), "COORD" => coord,
                             "ZCORN" => zcorn, "ACTNUM" => fill(true, nx, ny, 1))
    @compile_workload begin
        # as atlas_mesh: the largest piece, meshed, its columns
        grid["ACTNUM"] = _largest_piece(grid)
        mesh = _mesh_from_grid_section(grid)
        columns, _ = ve_columns(mesh)
        s = trap_surface(mesh, columns)
        for method in (:node, :cell)
            res = trap_analysis(s; method)
            trap_volumes(s, res)
            maximum(res.traps) > 0 && downstream_traps(res, 1)
        end
    end
end

# Comparisons of the CO2 Storage Atlas code with co2lab (roadmap Stage 5), used
# by the package tests (runtests.jl) and the gate (scripts/stage5_gate.jl).
# Each takes an MRST export (scripts/mrst_dump_atlas_*.m, read with matread)
# and returns an OracleReport.

# atlas_maps.mat (scripts/mrst_dump_atlas_maps.m): every Atlas map as co2lab's
# readAAIGrid reads it, getAtlasGrid's maps for five formations at coarsening 1
# and 3, and their average rock. Julia reads with co2lab's row handling here;
# its default reading is the same map one row further south.
function _atlas_vs_mrst(d)
    rep = OracleReport("Atlas maps vs co2lab"; mrst = d)
    # NaN count, sum and position-weighted sum, as the dump's fingerprint
    weighted(x) = sum(x[i, j] * (i + 7j) for i in axes(x, 1), j in axes(x, 2)
                      if !isnan(x[i, j]))
    fp(x) = (count(isnan, x), sum(filter(!isnan, x)), weighted(x))
    exact = Tolerance(0.0, 0.0)

    for e in d["files"]
        path = joinpath(atlas_data_directory(), e["dataset"], e["file"])
        g = read_aai_grid(path; co2lab_rows = true)
        m = g.meta
        label = e["file"]
        check!(rep, "$label: header", [m.ncols, m.nrows, m.xllcorner, m.yllcorner,
                                      m.cellsize, m.nodata], vec(e["meta"]); tol = exact)
        f = fp(g.data)
        check!(rep, "$label: NaN count", f[1], e["nan_count"]; tol = exact)
        # summation order differs (MATLAB's sum, Julia's pairwise sum)
        check!(rep, "$label: sums", [f[2], f[3]], [e["sum"], e["weighted_sum"]];
               tol = Tolerance(1e-12, 0.0))
        if !isempty(e["data"])
            same = isequal(g.data, e["data"])
            check!(rep, "$label: values", same ? 0 : 1, 0; tol = exact)
        end
    end

    for (key, maps) in d["maps"]
        name, c = split(key, "_c")
        ours = read_atlas_maps([name, "Jurassic", "Jurassicmid"];
                               coarsening = parse(Int, c), co2lab_compatible = true)
        for e in maps
            o = only(filter(x -> x.name == e["name"] &&
                                 string(x.variant) == e["variant"], ours))
            f = fp(o.data)
            check!(rep, "$key $(e["name"]) $(e["variant"])",
                   [f[1], f[2], f[3], o.meta.cellsize, o.meta.dims...],
                   [e["nan_count"], e["sum"], e["weighted_sum"], e["cellsize"],
                    vec(e["dims"])...]; tol = Tolerance(1e-12, 0.0))
        end
    end

    # getAtlasGrid replaces the table's values by the averages of heterogeneous
    # rock maps where a formation has them (Sto): step 2
    rock_suffixes = ("perm", "poro", "ntg")
    with_rock_maps = Set(f.name[1:end-length(x)] for f in atlas_map_files()
                         for x in rock_suffixes if endswith(f.name, x))
    for (name, v) in d["petro"]
        name in with_rock_maps && continue
        r = atlas_average_rock(name)
        check!(rep, "$name: average perm, poro", [r.perm, r.poro], vec(v); tol = exact)
    end
    return rep
end

# atlas_grids.mat (scripts/mrst_dump_atlas_grids.m): getAtlasGrid's decks and
# petroinfo for ten formation/coarsening/refining/layer cases. Julia builds
# them with co2lab_compatible = true; co2lab's Inf at undefined pillars is 0
# here (GeoEnergyIO's convention).
function _atlas_grids_vs_mrst(d)
    rep = OracleReport("Atlas grid sections vs co2lab"; mrst = d)
    exact = Tolerance(0.0, 0.0)
    same(a, b) = isequal(a, b) ? 0 : 1
    for e in d["cases"]
        c, r, nz = Int(e["coarsening"]), Int(e["refining"]), Int(e["nz"])
        o = atlas_grid_section(e["name"]; coarsening = c, refining = r, nz,
                               min_thickness = 0.0, co2lab_compatible = true)
        g = o.grid
        label = "$(e["name"]) c$c r$r nz$nz"
        zero_inf(x) = (v = vec(x); ifelse.(isinf.(v), 0.0, v))

        check!(rep, "$label: cartDims", collect(g["cartDims"]), vec(e["cartDims"]);
               tol = exact)
        for k in ("COORD", "ZCORN")
            check!(rep, "$label: $k", same(g[k], zero_inf(e[k])), 0; tol = exact)
        end
        check!(rep, "$label: ACTNUM", same(vec(g["ACTNUM"]), vec(e["ACTNUM"] .> 0)), 0;
               tol = exact)
        # co2lab's rock arrays have one layer (co2lab-issues.org)
        for k in ("PORO", "PERMX", "NTG")
            haskey(e, k) || continue
            check!(rep, "$label: $k", same(vec(g[k][:, :, 1]), vec(e[k])), 0;
                   tol = exact)
        end
        # petroinfo [avgperm, avgporo, avgntg], avgntg NaN where co2lab has none
        p, q = o.petro, vec(e["petro"])
        ours = [p.perm, p.poro, something(p.ntg, NaN)]
        has = .!isnan.(q)
        check!(rep, "$label: petroinfo defined", same(.!isnan.(ours), has), 0;
               tol = exact)
        check!(rep, "$label: petroinfo", ours[has], q[has]; tol = Tolerance(1e-12, 0.0))
    end
    return rep
end

# atlas_formations.mat (scripts/mrst_dump_atlas_formations.m):
# getFormationTopGrid's Gt and rock2D for eight formation/coarsening/layer
# cases. Julia removes cells thinner than 1 mm, which co2lab keeps: its
# columns must be co2lab's minus thin ones, matched by the logical index of
# their top cell.
function _atlas_domains_vs_mrst(d)
    rep = OracleReport("Atlas VE domains vs co2lab (getFormationTopGrid)"; mrst = d)
    exact = Tolerance(0.0, 0.0)
    for e in d["cases"]
        c, L = Int(e["coarsening"]), Int(e["layers"])
        o = atlas_ve_domain(e["name"]; coarsening = c, layers = L,
                            co2lab_compatible = true)
        m, cols, ve = o.mesh, o.columns, o.ve
        label = "$(e["name"]) c$c layers $L"
        LI = LinearIndices(Tuple(Int.(vec(e["cartDims"]))))
        lin(cell) = LI[cell_ijk(m, cell)...]

        mrst_cells = Set(Int.(vec(e["indexMap"])))
        check!(rep, "$label: 3D cells not in co2lab's",
               count(i -> !(lin(i) in mrst_cells), 1:number_of_cells(m)), 0; tol = exact)
        top = Dict(t => k for (k, t) in enumerate(Int.(vec(e["top_cell"]))))
        p = [get(top, lin(cols.cells[cols.pos[k]]), 0)
             for k in 1:number_of_columns(cols)]
        check!(rep, "$label: columns not in co2lab's", count(==(0), p), 0; tol = exact)
        any(==(0), p) && continue

        # co2lab's columns that are not here: all thinner than 1 mm per layer
        H_m = vec(e["height"])
        gone = setdiff(1:Int(e["num_cells"]), p)
        check!(rep, "$label: co2lab's columns removed (recorded)", length(gone),
               length(gone); tol = exact)
        check!(rep, "$label: removed columns thinner than 1 mm per layer",
               count(>=(1e-3 * L), H_m[gone]), 0; tol = exact)

        dom, geo, rock = ve.domain, ve.geometry, ve.rock
        N = dom[:neighbors]
        F = Int.(e["face_pairs"])
        kept = Set(p)
        pairs_m = Set(minmax(F[f, 1], F[f, 2]) for f in axes(F, 1)
                      if F[f, 1] in kept && F[f, 2] in kept)
        same_faces = Set(minmax(p[N[1, f]], p[N[2, f]]) for f in axes(N, 2)) == pairs_m
        check!(rep, "$label: interior faces", same_faces ? 0 : 1, 0; tol = exact)

        rel = Tolerance(1e-12, 0.0)
        check!(rep, "$label: caprock depth", geo.z_top, vec(e["z_caprock"])[p];
               tol = rel)
        check!(rep, "$label: height", geo.H, H_m[p]; tol = Tolerance(1e-10, 1e-9))
        check!(rep, "$label: perm", rock.perm, vec(e["perm"])[p]; tol = rel)
        check!(rep, "$label: poro", rock.poro, vec(e["poro"])[p]; tol = rel)
        has_ntg = !isempty(e["ntg"])
        check!(rep, "$label: has ntg", isnothing(rock.ntg) ? 0 : 1, has_ntg ? 1 : 0;
               tol = exact)
        has_ntg && check!(rep, "$label: ntg", rock.ntg, vec(e["ntg"])[p]; tol = rel)

        # pore volume: 3D cell volumes against MRST's area·H (Stage 1, ≤ 0.2 %)
        ntg = has_ntg ? rock.ntg : 1.0
        pv = dom[:volumes] .* dom[:porosity] .* ntg
        check!(rep, "$label: pore volume (≤ 0.2 %)", pv, vec(e["pore_volume"])[p];
               tol = Tolerance(2e-3, 0.0))
    end
    return rep
end

# atlas_all.mat (scripts/mrst_dump_atlas_all.m): the formations getAtlasGrid
# builds decks for, and getFormationTopGrid's column count, total height,
# total pore volume and mean rock for each at coarsening 5. Julia removes
# columns thinner than 1 mm; formations without Atlas rock get placeholder
# values (their rock is not compared).
function _atlas_all_vs_mrst(d)
    rep = OracleReport("All Atlas formations vs co2lab (coarsening 5)"; mrst = d)
    exact = Tolerance(0.0, 0.0)
    check!(rep, "formations with a grid", Set(atlas_grid_formations()) ==
           Set(d["deck_names"]) ? 0 : 1, 0; tol = exact)
    for e in d["cases"]
        name = e["name"]
        norock = isnan(e["mean_perm"]) || isnan(e["mean_poro"])
        kw = norock ? (perm = 1e-13, poro = 0.2) : (;)
        o = atlas_ve_domain(name; coarsening = Int(e["coarsening"]),
                            co2lab_compatible = true, kw...)
        ve, n = o.ve, number_of_columns(o.columns)

        check!(rep, "$name: columns removed (≤ 1.5 %)", (Int(e["num_cells"]) - n) / n,
               0.0; tol = Tolerance(0.0, 0.015))
        check!(rep, "$name: total height", sum(ve.geometry.H), e["sum_height"];
               tol = Tolerance(1e-7, 0.0))
        norock && continue
        ntg = isnothing(ve.rock.ntg) ? 1.0 : ve.rock.ntg
        pv = sum(ve.domain[:volumes] .* ve.domain[:porosity] .* ntg)
        check!(rep, "$name: total pore volume", pv, e["pore_volume"];
               tol = Tolerance(1e-4, 0.0))
        check!(rep, "$name: mean perm, poro", [mean(ve.rock.perm), mean(ve.rock.poro)],
               [e["mean_perm"], e["mean_poro"]]; tol = Tolerance(1e-12, 0.0))
    end
    return rep
end

# atlas_volumes.mat (scripts/mrst_dump_atlas_volumes.m): modelsFromAtlas.m's
# table for the North Sea formations at full resolution — columns, caprock
# depth range, bulk volume — as examples/atlas/models_from_atlas.jl prints it.
# co2lab keeps processGRDECL's first piece, Julia the largest; Julia removes
# cells thinner than 1 mm. Only in the gate (about a minute).
function _atlas_volumes_vs_mrst(d)
    rep = OracleReport("modelsFromAtlas table vs co2lab (full resolution)"; mrst = d)
    exact = Tolerance(0.0, 0.0)
    for e in d["cases"]
        name = e["name"]
        pieces = e["piece_cells"] isa Number ? [e["piece_cells"]] : vec(e["piece_cells"])
        check!(rep, "$name: co2lab's first piece is its largest",
               pieces[1] == maximum(pieces) ? 0 : 1, 0; tol = exact)

        m = atlas_mesh(name; co2lab_compatible = true)
        n = number_of_columns(m.columns)
        z_top = ve_column_geometry(m.mesh, m.columns).z_top
        volume = sum(tpfv_geometry(m.mesh).volumes)
        check!(rep, "$name: columns removed (≤ 1.5 %)", (e["cells"] - n) / n, 0.0;
               tol = Tolerance(0.0, 0.015))
        # caprock depth range: co2lab's when no column is removed; otherwise
        # within it (removing sub-millimetre cells can detach a thin-walled
        # lobe, e.g. Utsira's deepest part, which then is not the largest piece)
        zmin, zmax = extrema(z_top)
        if n == e["cells"]
            check!(rep, "$name: caprock depth range", [zmin, zmax],
                   [e["zmin"], e["zmax"]]; tol = Tolerance(1e-12, 0.0))
        else
            check!(rep, "$name: caprock depth range within co2lab's",
                   [max(e["zmin"] - zmin, 0.0), max(zmax - e["zmax"], 0.0)], [0.0, 0.0];
                   tol = Tolerance(0.0, 1e-9))
            check!(rep, "$name: depth range difference (recorded, m)",
                   [zmin - e["zmin"], zmax - e["zmax"]], [0.0, 0.0];
                   tol = Tolerance(0.0, Inf))
        end
        check!(rep, "$name: bulk volume", volume, e["volume"];
               tol = Tolerance(1e-4, 0.0))
    end
    return rep
end

# sleipner.mat (scripts/mrst_dump_sleipner.m): makeSleipnerVEmodel's steps with
# global coordinates — cells, columns, faces, caprock depths, heights, rock,
# pore volumes and the open boundary. 43 layers: only in the gate.
function _sleipner_vs_mrst(d)
    rep = OracleReport("Sleipner VE model vs co2lab (makeSleipnerVEmodel)"; mrst = d)
    exact = Tolerance(0.0, 0.0)
    s = sleipner_ve_domain()
    m, cols, ve = s.mesh, s.columns, s.ve
    LI = LinearIndices(Tuple(Int.(vec(d["cartDims"]))))
    lin(c) = LI[cell_ijk(m, c)...]

    cells = sort([lin(c) for c in 1:number_of_cells(m)])
    check!(rep, "3D cells", cells == sort(Int.(vec(d["indexMap"]))) ? 0 : 1, 0;
           tol = exact)
    top = Dict(t => k for (k, t) in enumerate(Int.(vec(d["top_cell"]))))
    p = [get(top, lin(cols.cells[cols.pos[k]]), 0) for k in 1:number_of_columns(cols)]
    check!(rep, "columns matched", count(>(0), p), Int(d["num_cells"]); tol = exact)
    any(==(0), p) && return rep

    N = ve.domain[:neighbors]
    F = Int.(d["face_pairs"])
    same = Set(minmax(p[N[1, f]], p[N[2, f]]) for f in axes(N, 2)) ==
           Set(minmax(F[f, 1], F[f, 2]) for f in axes(F, 1))
    check!(rep, "interior faces", same ? 0 : 1, 0; tol = exact)
    check!(rep, "open boundary sides (bcIxVE)", length(s.sides),
           Int(d["num_boundary_faces"]); tol = exact)

    geo, rock = ve.geometry, ve.rock
    rel = Tolerance(1e-12, 0.0)
    check!(rep, "caprock depth", geo.z_top, vec(d["z_caprock"])[p]; tol = rel)
    check!(rep, "height", geo.H, vec(d["height"])[p]; tol = Tolerance(1e-10, 1e-9))
    check!(rep, "perm", rock.perm, vec(d["perm"])[p]; tol = rel)
    check!(rep, "poro", rock.poro, vec(d["poro"])[p]; tol = rel)

    # global coordinates: the top cells' centroids
    g = tpfv_geometry(m)
    tops = [cols.cells[cols.pos[k]] for k in 1:number_of_columns(cols)]
    check!(rep, "top cell centroids (global x, y, z)", vec(g.cell_centroids[:, tops]'),
           vec(d["top_centroids"][p, :]); tol = rel)
    # column volumes from the 3D cells, on both sides (with porosity checked
    # above, this fixes the pore volume)
    check!(rep, "column volumes (3D cells)", ve.domain[:volumes],
           vec(d["column_volume"])[p]; tol = Tolerance(1e-10, 0.0))
    # pore volume = column volume x thickness-weighted porosity (co2lab's
    # averageRock); recorded: against co2lab's VE pore volume (top-face area
    # x H, by design) and against the volume-weighted sum over the 3D cells
    pv = ve.domain[:volumes] .* ve.domain[:porosity]
    check!(rep, "co2lab's area x H pore volume (recorded)", pv,
           vec(d["pore_volume"])[p]; tol = Tolerance(Inf, 0.0))
    check!(rep, "sum of 3D cell pore volumes (recorded)", pv,
           vec(d["pore_volume_3d"])[p]; tol = Tolerance(Inf, 0.0))
    check!(rep, "bulk volume", sum(g.volumes), d["bulk_volume"]; tol = rel)
    return rep
end

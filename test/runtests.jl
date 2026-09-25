using Test
using JutulCO2Lab
using JutulCO2Lab.Oracle
using Jutul
using JutulDarcy
using MAT

@testset "JutulCO2Lab" begin

    @testset "Oracle" begin
        tol = Oracle.Tolerance(1e-10, 1e-12)

        @testset "identical input passes" begin
            rep = OracleReport("identical")
            a = [1.0, 2.0, 3.0]
            check!(rep, "a", a, copy(a); tol = tol)
            @test passed(rep)
        end

        @testset "a difference above tolerance fails" begin
            rep = OracleReport("perturbed")
            check!(rep, "a", [1.0, 2.0, 3.0], [1.0, 2.0, 3.1]; tol = tol)
            @test !passed(rep)
            @test rep.results[1].maxabs ≈ 0.1
        end

        @testset "an empty report is not a pass" begin
            @test !passed(OracleReport("nothing ran"))
        end

        @testset "an empty comparison is not a pass" begin
            rep = OracleReport("empty")
            check!(rep, "a", Float64[], Float64[]; tol = tol)
            @test !passed(rep)
        end

        @testset "length mismatch fails without throwing" begin
            rep = OracleReport("ragged")
            check!(rep, "a", [1.0, 2.0], [1.0, 2.0, 3.0]; tol = tol)
            @test !passed(rep)
            @test occursin("length mismatch", rep.results[1].note)
        end

        @testset "non-finite entries fail" begin
            rep = OracleReport("nan")
            check!(rep, "a", [1.0, NaN], [1.0, 2.0]; tol = tol)
            @test !passed(rep)
            @test occursin("non-finite", rep.results[1].note)
        end

        @testset "perm reorders the reference" begin
            rep = OracleReport("permuted")
            check!(rep, "a", [3.0, 1.0, 2.0], [1.0, 2.0, 3.0];
                   tol = tol, perm = [3, 1, 2])
            @test passed(rep)
        end

        @testset "MATLAB shapes flatten" begin
            rep = OracleReport("shapes")
            check!(rep, "scalar", 2.0, fill(2.0, 1, 1); tol = tol)   # MATLAB 1x1
            check!(rep, "column", [1.0, 2.0], reshape([1.0, 2.0], 2, 1); tol = tol)
            @test passed(rep)
        end

        @testset "provenance" begin
            srcs = source_provenance()
            @test [s.name for s in srcs] == ["Jutul", "JutulDarcy"]
            for s in srcs
                @test isdir(s.dir)
                @test s.commit != "not loaded"
                # Both are git checkouts here; an 8-char hex hash, not "unknown".
                @test occursin(r"^[0-9a-f]{8}$", s.commit)
                @test s.dirty isa Bool
            end

            # A report records them, plus MRST when a dump is supplied.
            fake_dump = Dict("mrst_commit" => "0123456789abcdef",
                             "mrst_root" => "/x/mrst")
            rep = OracleReport("with mrst"; mrst = fake_dump)
            @test length(rep.sources) == 3
            @test rep.sources[end].name == "MRST"
            @test rep.sources[end].commit == "01234567"
            @test occursin("Jutul ", sprint(show, MIME"text/plain"(), rep))

            # A dump predating provenance still works.
            @test Oracle._mrst_source_info(Dict{String, Any}()).commit == "unknown"

            # Not a git repo: reported, not thrown.
            info = Oracle._git_source_info("tmp", mktempdir())
            @test info.commit == "unknown"
            @test info.dirty === missing
        end

        @testset "Stage 1 tolerances are defined for every gate quantity" begin
            for k in (:height, :pore_volume, :transmissibility, :trans_mult, :well_index)
                @test haskey(STAGE1_TOL, k)
            end
        end
    end

    # VE columns and column geometry against MRST reference data. The .mat
    # files are generated (scripts/run_mrst_dump_octave.sh, mrst_dump_grids.m,
    # mrst_dump_holes_case.m) and not versioned, so these skip if absent.
    refdir = joinpath(@__DIR__, "..", "..", "analysis", "reference-data")
    iv(x) = x isa Number ? [Int(x)] : Int.(vec(x))
    exact = Tolerance(0.0, 0.0)
    geomtol = Tolerance(1e-12, 1e-9)

    @testset "VE columns: Johansen vs MRST topSurfaceGrid" begin
        f = joinpath(refdir, "johansen_grids.mat")
        if !isfile(f)
            @info "Skipping: $f not generated"
        else
            d = matread(f); Gt = d["Gt"]
            mesh = ve_mesh_from_mrst(d["G3"])
            cols, discarded = ve_columns(mesh)
            geo = ve_column_geometry(mesh, cols)
            rep = OracleReport("VE columns, Johansen")
            check!(rep, "discarded", length(discarded), 0; tol = exact)
            check!(rep, "columns.cells", cols.cells, iv(Gt["columns"]["cells"]);
                   tol = exact)
            check!(rep, "columnPos", cols.pos, iv(Gt["cells"]["columnPos"]); tol = exact)
            check!(rep, "caprock z", geo.z_top, Gt["cells"]["z"]; tol = geomtol)
            check!(rep, "H", geo.H, Gt["cells"]["H"]; tol = geomtol)
            check!(rep, "dz", geo.dz, Gt["columns"]["dz"]; tol = geomtol)
            check!(rep, "z", geo.z, Gt["columns"]["z"]; tol = geomtol)
            passed(rep) || display(rep)
            @test passed(rep)
        end
    end

    @testset "VE columns: discard_below_holes vs MRST" begin
        f = joinpath(refdir, "holes_case.mat")
        if !isfile(f)
            @info "Skipping: $f not generated"
        else
            d = matread(f)
            mesh = ve_mesh_from_mrst(d["G_in"])
            rep = OracleReport("VE columns, holes case")
            ck, dk = ve_columns(mesh; discard_below_holes = false)
            check!(rep, "discard=false: discarded", length(dk), 0; tol = exact)
            check!(rep, "discard=false: columns.cells", ck.cells,
                   iv(d["keep"]["columns_cells"]); tol = exact)
            _, disc = ve_columns(mesh)
            # MRST lists discarded cells in column-major matrix order; compare as sets
            check!(rep, "discarded cells", disc, sort(iv(d["discard"]["discarded"]));
                   tol = exact)
            trimmed = ve_trim(mesh, disc)
            check!(rep, "trimmed logical index map", trimmed.cell_map,
                   iv(d["discard"]["indexMap_out"]); tol = exact)
            c2, disc2 = ve_columns(trimmed)
            check!(rep, "no further discards", length(disc2), 0; tol = exact)
            check!(rep, "columns.cells", c2.cells, iv(d["discard"]["columns_cells"]);
                   tol = exact)
            check!(rep, "columnPos", c2.pos, iv(d["discard"]["columnPos"]); tol = exact)
            g = ve_column_geometry(trimmed, c2)
            check!(rep, "H", g.H, d["discard"]["H"]; tol = geomtol)
            check!(rep, "caprock z", g.z_top, d["discard"]["z"]; tol = geomtol)
            passed(rep) || display(rep)
            @test passed(rep)
            coarse = CoarseMesh(trimmed, c2.partition)
            @test number_of_cells(coarse) == number_of_columns(c2)
        end
    end

    @testset "VE domain and boundary sides: Johansen vs MRST" begin
        fg = joinpath(refdir, "johansen_grids.mat")
        fs = joinpath(refdir, "stage1_johansen.mat")
        if !(isfile(fg) && isfile(fs))
            @info "Skipping: Johansen reference data not generated"
        else
            d = matread(fg); s1 = matread(fs)
            mesh = ve_mesh_from_mrst(d["G3"])
            cols, _ = ve_columns(mesh)
            r3 = d["rock3D"]
            K = Float64.(r3["perm"])
            ve = @test_logs ve_domain(mesh, cols; perm_x = K[:, 1], perm_y = K[:, 2],
                                      poro = vec(Float64.(r3["poro"])))
            sides, _ = ve_boundary_sides(mesh, cols, ve.domain)
            rep = OracleReport("VE domain, Johansen")
            check!(rep, "perm", ve.rock.perm, s1["perm"]; tol = Tolerance(1e-12, 0.0))
            check!(rep, "poro", ve.rock.poro, s1["poro"]; tol = Tolerance(1e-12, 0.0))
            pv = ve.domain[:volumes] .* ve.domain[:porosity]
            check!(rep, "pore volume", pv, s1["pore_volume"]; tol = Tolerance(2e-3, 0.0))
            zt = vec(s1["z_caprock"]); N = ve.domain[:neighbors]
            gdz_mrst = [-Jutul.gravity_constant * (zt[N[2, f]] - zt[N[1, f]])
                        for f in 1:size(N, 2)]
            check!(rep, "gdz", ve.domain[:gdz], gdz_mrst;
                   tol = Tolerance(1e-12, 1e-9))
            Nt = Int.(d["Gt"]["faces"]["neighbors"])
            n_bnd_edges = count((Nt[:, 1] .== 0) .| (Nt[:, 2] .== 0))
            check!(rep, "boundary sides", length(sides), n_bnd_edges; tol = exact)
            passed(rep) || display(rep)
            @test passed(rep)
        end
    end

    @testset "face directions and column sides: Johansen (topology)" begin
        f = joinpath(refdir, "johansen_grids.mat")
        if !isfile(f)
            @info "Skipping: $f not generated"
        else
            d = matread(f); G = d["G3"]
            mesh = ve_mesh_from_mrst(G)
            cols, _ = ve_columns(mesh)
            # the logical-index rule reproduces MRST's tag for every interior cell-face
            N = Int.(G["faces"]["neighbors"]); cf = G["cells"]["faces"]
            fp = Int.(vec(G["cells"]["facePos"]))
            ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
            # cell-face row p of cell c: missing if boundary, else does the rule agree?
            function tag_agrees(c, p)
                f = Int(cf[p, 1]); l, r = N[f, 1], N[f, 2]
                (l > 0 && r > 0) || return missing
                o = l == c ? r : l
                return JutulCO2Lab._logical_direction(ijk[c], ijk[o]) == cf[p, 2]
            end
            agree = skipmissing(tag_agrees(c, p) for c in 1:(length(fp) - 1)
                                for p in fp[c]:(fp[c+1] - 1))
            agree = collect(agree)
            @test length(agree) == 159578 && all(agree)
            # column sides: neighbours are exactly the CoarseMesh's column pairs
            r3 = d["rock3D"]
            K = Float64.(r3["perm"])
            ve = @test_logs ve_domain(mesh, cols; perm_x = K[:, 1], perm_y = K[:, 2],
                                      poro = vec(Float64.(r3["poro"])))
            cs = ve_column_sides(mesh, cols)
            Nc = ve.domain[:neighbors]
            side_pairs = Set(minmax(cs.column[i], cs.neighbour[i])
                             for i in 1:length(cs) if cs.neighbour[i] != 0)
            @test side_pairs == Set(minmax(Nc[1, f], Nc[2, f]) for f in 1:size(Nc, 2))
            # as many outer sides as MRST has boundary edges
            @test count(==(0), cs.neighbour) == 456
            sides, unassigned = ve_boundary_sides(mesh, cols, ve.domain)
            # exposed faces of internal steps
            @test length(unassigned) == 718
            # fully offset faults: logical neighbour exists, no connection
            @test count(sides.across_fault) == 88
            fault_i = Set(cols.pillar[c][1] for c in sides.column[sides.across_fault])
            @test fault_i == Set([43, 44, 78, 79])
        end
        @testset "Jutul-built corner-point grid with a fault (no MRST)" begin
            # 3x2x3 cells, 100 m x 10 m; the i = 3 columns shifted down 15 m
            nx, ny, nz = 3, 2, 3
            coord = Float64[]
            for j in 0:ny, i in 0:nx
                append!(coord, [i * 100.0, j * 100.0, 0.0, i * 100.0, j * 100.0, 2000.0])
            end
            zc = [1000.0 + (k - 1 + s) * 10.0 + ((ii + 1) ÷ 2 == 3 ? 15.0 : 0.0)
                  for k in 1:nz for s in (0, 1) for jj in 1:2ny for ii in 1:2nx]
            grid = Dict("cartDims" => [nx, ny, nz], "COORD" => coord, "ZCORN" => zc,
                        "ACTNUM" => fill(true, nx, ny, nz))
            mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
            cols, disc = ve_columns(mesh)
            @test number_of_columns(cols) == 6 && isempty(disc)
            geo = ve_column_geometry(mesh, cols)
            @test geo.z_top ≈ [1000, 1000, 1015, 1000, 1000, 1015] && all(geo.H .≈ 30)
            cs = ve_column_sides(mesh, cols)
            # outer sides: the perimeter
            @test count(==(0), cs.neighbour) == 2 * (nx + ny)
            side_pairs = Set(minmax(cs.column[i], cs.neighbour[i])
                             for i in 1:length(cs) if cs.neighbour[i] != 0)
            @test length(side_pairs) == (nx - 1) * ny + nx * (ny - 1)
            ve = ve_domain(mesh, cols; perm_x = fill(1e-13, 18), poro = fill(0.2, 18))
            # lateral anisotropy: warning, and the VE permeability is x's
            kx, ky = fill(1e-13, 18), fill(2e-13, 18)
            rock = @test_logs (:warn, r"anisotropic in 18 of 18") ve_average_rock(
                cols, geo.dz; perm_x = kx, perm_y = ky, poro = fill(0.2, 18))
            @test rock.perm ≈ fill(1e-13, 6)
            sides, unassigned = ve_boundary_sides(mesh, cols, ve.domain)
            @test length(sides) == 10
            @test sort(unique(round.(sides.edge_depth, digits = 6))) == [1000.0, 1015.0]
            # partial offset: still connected, so no sealed-fault sides
            @test !any(sides.across_fault)
            # The fault's exposed strip above i = 3: upper part of the i = 2 columns'
            # right side (cells k = 1, 2), tagged :right by Jutul, no-flow. (The strip
            # below i = 2 on the i = 3 columns' left side is not represented as faces
            # by GeoEnergyIO at all — MRST's grids have it; cell volumes are
            # unaffected.)
            bdir = JutulCO2Lab._boundary_directions(mesh)
            @test length(unassigned) == 4 && all(==(2), bdir[unassigned])
            @test all(cell_ijk(mesh, mesh.boundary_faces.neighbors[f])[1] == 2
                      for f in unassigned)
        end
        dir = JutulCO2Lab._logical_direction
        @test_throws ErrorException dir((1, 1, 1), (2, 2, 1))   # diagonal
        @test dir((3, 4, 5), (3, 4, 2)) == 5                    # top
        @test dir((3, 4, 5), (4, 4, 9)) == 2                    # i+, across a fault
    end

    # Proves the toolchain end-to-end: local Jutul + JutulDarcy assemble,
    # linearise and solve. Stock immiscible physics, no VE.
    @testset "toolchain smoke test" begin
        Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)

        mesh = reservoir_mesh((4, 4, 2), (200.0, 200.0, 20.0))
        domain = reservoir_domain(mesh, permeability = 0.1Darcy, porosity = 0.2)
        inj = setup_vertical_well(domain, 1, 1, name = :Injector)
        prod = setup_vertical_well(domain, 4, 4, name = :Producer)

        rhoLS, rhoGS = 1000.0kg/meter^3, 100.0kg/meter^3
        sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()),
                               reference_densities = [rhoLS, rhoGS])
        model, parameters = setup_reservoir_model(domain, sys,
                                                  wells = [inj, prod], extra_out = true)
        state0 = setup_reservoir_state(model, Pressure = 150bar,
                                       Saturations = [1.0, 0.0])

        dt = fill(30.0day, 5)
        pv = pore_volume(model, parameters)
        rate = 0.1 * sum(pv) / sum(dt)
        controls = Dict(
            :Injector => InjectorControl(TotalRateTarget(rate), [0.0, 1.0],
                                         density = rhoGS),
            :Producer => ProducerControl(BottomHolePressureTarget(100bar)),
        )
        forces = setup_reservoir_forces(model, control = controls)

        result = simulate_reservoir(state0, model, dt,
                                    parameters = parameters, forces = forces,
                                    info_level = -1)

        @test length(result.states) == length(dt)
        sg = result.states[end][:Saturations][2, :]
        @test maximum(sg) > 0.0                 # gas entered the reservoir
        @test all(isfinite, result.states[end][:Pressure])
    end


    @testset "co2lab fluid system (Stage 2 item 0)" begin
        Darcy, day = si_units(:darcy, :day)
        # 1D sloping aquifer, injection then shut-in (the prototype's case, shortened)
        nx = 20
        mesh = UnstructuredMesh(CartesianMesh((nx, 1, 1), (1000.0, 100.0, 20.0)))
        for (k, x) in enumerate(mesh.node_points)
            mesh.node_points[k] = x .+ (0.0, 0.0, 1000.0 + 0.02 * x[1])
        end
        domain = reservoir_domain(mesh, permeability = 0.5Darcy, porosity = 0.25)
        inj = setup_well(domain, [nx - 3]; name = :inj)
        p0 = [1020.0 * 9.81 * c[3] for c in eachcol(domain[:cell_centroids])]
        q_in = 1e4 / day
        n_on, n_off = 10, 10
        dt = fill(10.0day, n_on + n_off)

        function run(; kw...)
            model = setup_co2lab_model(domain; wells = [inj], kw...)
            state0 = setup_co2lab_state(model; pressure = p0)
            on = InjectorControl(TotalRateTarget(q_in), [0.0, 1.0], density = 1.87)
            f_on = setup_reservoir_forces(model, control = Dict(:inj => on))
            off = DisabledControl()
            f_off = setup_reservoir_forces(model, control = Dict(:inj => off))
            forces = vcat(fill(f_on, n_on), fill(f_off, n_off))
            case = JutulCase(model, dt, forces; state0 = state0,
                             parameters = setup_parameters(model))
            sim, cfg = setup_co2lab_simulator(case; info_level = -1)
            return model, simulate_reservoir(case; simulator = sim, config = cfg).states
        end
        injected = n_on * 10day * q_in * 1.87
        co2_mass(s) = sum(s[:TotalMasses][2, :])

        @testset "configuration errors" begin
            # hysteresis: Stage 3; dissolution needs a rate; instantaneous: Stage 5
            @test_throws ErrorException setup_co2lab_model(domain; hysteresis = true)
            @test_throws ErrorException setup_co2lab_model(domain; rs_sat = 25.0)
            @test_throws ErrorException setup_co2lab_model(domain; rs_sat = 25.0,
                                                           transfer_rate = Inf)
        end

        @testset "sGmax update: running maximum, capped at 1 - res_water" begin
            v = MaxGasSaturation(0.2)
            sGmax = [0.0, 0.5, 0.7]
            JutulCO2Lab._update_max_gas_saturation!(sGmax, v, [0.3, 0.4, 0.9])
            @test sGmax == [0.3, 0.5, 0.8]
        end

        @testset "dissolution off: unknowns p, sG" begin
            model, st = run()
            r = model[:Reservoir]
            @test collect(keys(r.primary_variables)) == [:Pressure, :sG]
            @test collect(keys(r.equations)) == [:mass_conservation]
            @test !has_dissolution(r.system)
            @test haskey(r.parameters, :sGmax)
            @test Jutul.matrix_layout(r.context) isa Jutul.BlockMajorLayout
            @test length(st) == length(dt)
            @test abs(co2_mass(st[end]) - injected) / injected < 1e-3
            runmax = accumulate((a, b) -> max.(a, b), [s[:sG] for s in st])
            @test all(s[:sGmax] ≈ m for (s, m) in zip(st, runmax))
        end

        @testset "dissolution on (rate-driven): unknowns p, sG, rs" begin
            rs_sat = 25.0
            model, st = run(rs_sat = rs_sat, transfer_rate = 1 / (60day))
            r = model[:Reservoir]
            @test collect(keys(r.primary_variables)) == [:Pressure, :sG, :Rs]
            @test haskey(r.equations, :dissolved_co2)
            @test has_dissolution(r.system)
            @test Jutul.matrix_layout(r.context) isa Jutul.BlockMajorLayout
            # wellbore at equilibrium: black-oil's own unknowns
            @test haskey(model[:inj].primary_variables, :BlackOilUnknown)
            @test length(st) == length(dt)
            @test abs(co2_mass(st[end]) - injected) / injected < 1e-3
            sG = reduce(hcat, [s[:sG] for s in st])
            Rs = reduce(hcat, [s[:Rs] for s in st])
            # free CO2 over undersaturated brine: impossible at equilibrium
            @test any((sG .> 1e-4) .& (Rs .< 0.99 * rs_sat))
            @test maximum(Rs) <= rs_sat * (1 + 1e-6)
            dis_frac(s) = sum(s[:DissolvedCO2]) / co2_mass(s)
            # dissolution continues after shut-in
            @test dis_frac(st[end]) > dis_frac(st[n_on])
        end
    end

end

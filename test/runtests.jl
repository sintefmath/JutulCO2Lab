using Test
using JutulCO2Lab
using JutulCO2Lab.Oracle
using Jutul
using JutulDarcy
using MAT

# --------------------- Helpers: small corner-point grids -----------------------

# An nx x 1 x nz corner-point grid built by Jutul, 100 m cells in x; layer k of
# pillar i spans depths top(i, k)..top(i, k + 1).
function _cp_mesh(nx, nz, top; actnum = fill(true, nx, 1, nz))
    coord = Float64[]
    for j in 0:1, i in 0:nx
        append!(coord, [i * 100.0, j * 100.0, 0.0, i * 100.0, j * 100.0, 2000.0])
    end

    zc = [top((ii + 1) ÷ 2, k + s) for k in 1:nz for s in (0, 1)
          for jj in 1:2 for ii in 1:2nx]
    grid = Dict("cartDims" => [nx, 1, nz], "COORD" => coord, "ZCORN" => zc,
                "ACTNUM" => actnum)
    return JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
end

# --------------- Helpers: comparisons with MRST dumps and runs ----------------

# The VE sloping-aquifer domain of an MRST dump (3D grid G3 and rock).
function _ve_case(d)
    mesh = ve_mesh_from_mrst(d["G3"])
    cols, _ = ve_columns(mesh)
    ve = ve_domain(mesh, cols; perm_x = vec(d["rock"]["perm"]),
                   poro = vec(d["rock"]["poro"]))

    sides, _ = ve_boundary_sides(mesh, cols, ve.domain)
    return ve, cols, sides
end

# mrst_ref_{rate,bhp,src}.mat (scripts/stage2_mrst_comparison.jl): co2lab's
# sloping aquifer run by MRST on Julia's geometry, same report steps, tight
# tolerances. With a source term instead of the well, Julia matches to solver
# precision (measured 5e-10); with the well, to the wellbore-storage
# difference (measured BHP 1e-5, rate 3e-4, pressure 5e-6, sG 1.3e-5).
function _ve_well_vs_mrst_runs(d, ref_rate, ref_bhp, ref_src)
    ve, cols, sides = _ve_case(d)
    g, rho_w = d["gravity"], d["rhow"]
    fluid = co2lab_ve_fluid(krmax = (1.0, 1.0), co2_rho_ref = d["rhoc"],
                            wat_rho_ref = rho_w, co2_mu_ref = d["muco2"],
                            wat_mu_ref = d["muw"])
    rc = (reference_pressure = d["p_ref"], compressibility = d["cf_rock"])

    bc = ve_boundary_conditions(ve, sides, ve_select_sides(sides, cols;
                                directions = (:right,)); pressure = z -> rho_w * g * z)

    p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = rho_w, g = g)
    col, dt = Int(d["well"]["cells"]), vec(d["dT"])
    on = vec(d["step_control"]) .== 1

    function run(target)
        if target isa Real
            # a CO2 source term at this rate instead of the well
            model = setup_co2lab_model(ve.domain; fluid = fluid,
                                       rock_compressibility = rc)
            src(q) = [SourceTerm(col, q; fractional_flow = [0.0, 1.0])]
            f_on = setup_reservoir_forces(model; bc = bc, sources = src(target))
            f_off = setup_reservoir_forces(model; bc = bc, sources = src(0.0))
        else
            w = setup_ve_well(ve, col; name = :inj, radius = d["well"]["r"])
            model = setup_co2lab_model(ve.domain; wells = [w], fluid = fluid,
                                       rock_compressibility = rc)
            forces(c) = setup_reservoir_forces(model; bc = bc, control = Dict(:inj => c))
            f_on = forces(co2_injector_control(target, fluid))
            f_off = forces(DisabledControl())
        end

        case = JutulCase(model, dt, [i ? f_on : f_off for i in on];
                         state0 = setup_co2lab_state(model; pressure = p0),
                         parameters = setup_parameters(model))

        sim, cfg = setup_co2lab_simulator(case; info_level = -1, timesteps = :none,
                                          initial_dt = maximum(dt), max_timestep = Inf,
                                          tol_cnv = 1e-9, tol_mb = 1e-12,
                                          max_nonlinear_iterations = 50)

        res = simulate_reservoir(case; simulator = sim, config = cfg)
        @test all(length(r[:ministeps]) == 1 for r in res.result.reports)
        return res
    end

    rep = OracleReport("well model: Julia vs MRST runs, sloping aquifer")
    field(res, k) = reduce(hcat, [s[k] for s in res.states])

    res = run(d["injection_rate"])
    check!(rep, "source term: pressure", field(res, :Pressure), ref_src["pressure"];
           tol = Tolerance(1e-8, 0.0))
    check!(rep, "source term: sG", field(res, :sG), ref_src["sG"];
           tol = Tolerance(0.0, 1e-8))

    for (name, ref, target) in
            (("rate", ref_rate, TotalRateTarget(d["injection_rate"])),
             ("bhp", ref_bhp, BottomHolePressureTarget(p0[col] + 5e5)))
        res = run(target)
        check!(rep, "$name: pressure", field(res, :Pressure), ref["pressure"];
               tol = Tolerance(2e-5, 0.0))
        check!(rep, "$name: sG", field(res, :sG), ref["sG"]; tol = Tolerance(0.0, 1e-4))

        if name == "rate"
            check!(rep, "rate: injector BHP", res.wells[:inj][:bhp][on],
                   vec(ref["bhp"])[on]; tol = Tolerance(1e-4, 0.0))
        else
            q = abs.(res.wells[:inj][:grat][on])
            check!(rep, "bhp: injection rate", q, vec(ref["qGs"])[on];
                   tol = Tolerance(1e-3, 0.0))
            check!(rep, "bhp: cumulative injection", sum(q .* dt[on]),
                   sum(vec(ref["qGs"])[on] .* dt[on]); tol = Tolerance(1e-4, 0.0))
        end
    end

    passed(rep) || display(rep)
    @test passed(rep)
end

# ve_reconstruction.mat: co2lab's VEstates23D on sampled states, three fluid
# cases; and, per column, reconstructed CO2 volume = VE CO2 volume.
function _ve_reconstruction_vs_mrst(d)
    ve, cols, _ = _ve_case(d)
    pv3 = vec(d["rock"]["poro"]) .* tpfv_geometry(ve.coarse.parent).volumes
    pv_ve = ve.domain[:volumes] .* ve.domain[:porosity]

    rep = OracleReport("VE reconstruction vs MRST co2lab")
    for (c, res) in enumerate(((0.0, 0.0), (0.2, 0.0), (0.2, 0.15)))
        fl = co2lab_ve_fluid(residual = res, co2_rho_ref = 700, wat_rho_ref = 1050,
                             wat_rho_pvt = (4.3e-10, 100e5))
        r = d["case$c"]
        for k in axes(r["p"], 2)
            st = Dict(:Pressure => r["p"][:, k], :sG => r["sg"][:, k],
                      :sGmax => r["sgmax"][:, k])

            o = ve_reconstruct(st, ve, cols, fl)
            check!(rep, "case $c, state $k: pressure", o.pressure, r["p3D"][:, k];
                   tol = Tolerance(1e-14, 0.0))
            check!(rep, "case $c, state $k: sG", o.sG, r["s3D"][:, k];
                   tol = Tolerance(0.0, 1e-12))
            check!(rep, "case $c, state $k: sGmax", o.sGmax, r["smax3D"][:, k];
                   tol = Tolerance(0.0, 1e-12))

            # without residual CO2 the plume holds all of the column's CO2
            res[2] == 0 || continue
            v3 = [sum(pv3[column_cells(cols, j)] .* o.sG[column_cells(cols, j)])
                  for j in 1:number_of_columns(cols)]
            check!(rep, "case $c, state $k: CO2 volume per column", v3,
                   pv_ve .* r["sg"][:, k]; tol = Tolerance(1e-11, 1e-9))
        end
    end

    passed(rep) || display(rep)
    @test passed(rep)
end

# ve_bc.mat: sides, boundary pressures and transmissibilities against co2lab's
# open boundary, and the boundary fluxes against co2lab's on sampled states.
function _ve_bc_vs_mrst(d)
    ve, cols, sides = _ve_case(d)
    ix = ve_select_sides(sides, cols)
    g = d["gravity"]

    # explicit map, our side -> MRST's boundary face, by caprock-edge midpoint
    fc = d["face_centroid"]
    dist(s, r) = hypot(sides.edge_midpoint[s][1] - fc[r, 1],
                       sides.edge_midpoint[s][2] - fc[r, 2])
    perm = [argmin([dist(s, r) for r in axes(fc, 1)]) for s in ix]
    @test length(ix) == size(fc, 1) && allunique(perm)
    @test maximum(dist(s, r) for (s, r) in zip(ix, perm)) < 1e-9

    own = ve_boundary_conditions(ve, sides, ix; pressure = z -> 1050 * g * z)
    bcs = ve_boundary_conditions(ve, sides, ix; pressure = z -> 1050 * g * z,
                                 trans = vec(d["T"])[perm])

    rep = OracleReport("VE boundary conditions vs MRST co2lab")
    exact = Tolerance(0.0, 0.0)
    check!(rep, "column", [b.cell for b in bcs], d["face_cell"]; perm = perm,
           tol = exact)
    check!(rep, "edge depth", sides.edge_depth[ix], d["face_z"]; perm = perm,
           tol = Tolerance(1e-12, 1e-9))
    check!(rep, "boundary pressure", [b.pressure for b in bcs], d["bc_value"];
           perm = perm, tol = Tolerance(1e-12, 1e-6))

    # our half transmissibility, from the 3D grid (Stage 1 geometry)
    check!(rep, "half transmissibility", [b.trans for b in own], d["T"];
           perm = perm, tol = Tolerance(5e-3, 0.0))

    r = d["samples"]
    for (ph, p, mob, rho, q) in ((1, "pW", "mobW", "rhoW", "qResW"),
                                  (2, "pG", "mobG", "rhoG", "qResG"))
        flux = [JutulCO2Lab._ve_bc_volume_flux(b, r[p][b.cell, s], r[mob][b.cell, s],
                    r[rho][b.cell, s], r["mobW"][b.cell, s] + r["mobG"][b.cell, s],
                    b.sat[ph], g) for b in bcs, s in axes(r["p"], 2)]
        check!(rep, "boundary flux, phase $ph", vec(flux), vec(r[q][perm, :]);
               tol = Tolerance(1e-11, 1e-20))
    end

    passed(rep) || display(rep)
    @test passed(rep)

    # the samples exercise inflow and outflow of brine and outflow of CO2
    @test any(<(0), r["qResW"]) && any(>(0), r["qResW"]) && any(<(0), r["qResG"])
end

# Model level: hydrostatic equilibrium stays at rest with every side open;
# with injection, CO2 in place + exited + in the wellbore = injected.
function _ve_bc_runs(d)
    ve, cols, sides = _ve_case(d)
    g, rhow, rhoc = d["gravity"], d["rhow"], d["rhoc"]
    fluid = co2lab_ve_fluid(krmax = (1.0, 1.0), co2_rho_ref = rhoc, wat_rho_ref = rhow,
                            co2_mu_ref = d["muco2"], wat_mu_ref = d["muw"])
    rc = (reference_pressure = d["p_ref"], compressibility = d["cf_rock"])

    p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = rhow, g = g)
    hydro(z) = rhow * g * z
    all_open = ve_boundary_conditions(ve, sides, ve_select_sides(sides, cols);
                                      pressure = hydro)

    dt = vec(d["dT"])[1:8]

    function run(bcs, wells, control)
        model = setup_co2lab_model(ve.domain; wells = wells, fluid = fluid,
                                   rock_compressibility = rc)
        push!(model[:Reservoir].output_variables, :RelativePermeabilities,
              :PhaseViscosities, :PhaseMassDensities)
        forces = setup_reservoir_forces(model; control = control, bc = bcs)
        case = JutulCase(model, dt, forces; parameters = setup_parameters(model),
                         state0 = setup_co2lab_state(model; pressure = p0))
        sim, cfg = setup_co2lab_simulator(case; info_level = -1,
                                          output_substates = true)
        return simulate_reservoir(case; simulator = sim, config = cfg).result
    end

    rest = run(all_open, [], nothing)
    @test maximum(abs.(rest.states[end][:Reservoir][:Pressure] .- p0)) < 1e-6
    @test maximum(rest.states[end][:Reservoir][:sG]) == 0.0

    # injection, open updip end only (co2lab's cross-sectional setup)
    right = ve_select_sides(sides, cols; directions = (:right,))
    bcs = ve_boundary_conditions(ve, sides, right; pressure = hydro)
    @test [b.cell for b in bcs] == [Int(d["bc"]["cell"])]

    w = setup_ve_well(ve, Int(d["well"]["cells"]); name = :inj)
    ctrl = co2_injector_control(TotalRateTarget(d["injection_rate"]), fluid)
    res = run(bcs, [w], Dict(:inj => ctrl))

    sub, dts, _ = Jutul.expand_to_ministeps(res)
    rate(s) = s[:Facility][:TotalSurfaceMassRate][1]
    injected = sum(rate(s) * t for (s, t) in zip(sub, dts))

    # exited mass, (brine, CO2), sub-step by sub-step
    exited = sum(vec(sum(ve_boundary_mass_fluxes(bcs, s[:Reservoir]), dims = 2)) * t
                 for (s, t) in zip(sub, dts))

    last = res.states[end]
    total = sum(last[:Reservoir][:TotalMasses][2, :]) + last[:inj][:TotalMasses][2, 1]

    # brine is pushed out; the CO2 has not reached the boundary yet
    @test exited[1] > 0
    @test (total + exited[2]) / injected ≈ 1 rtol = 1e-5

    # dissolution: Stage 5
    model = setup_co2lab_model(ve.domain; fluid = fluid, rock_compressibility = rc,
                               rs_sat = 0.03,
                               transfer_rate = 1e-7)
    forces = setup_reservoir_forces(model; bc = bcs)
    case = JutulCase(model, dt[1:1], forces; parameters = setup_parameters(model),
                     state0 = setup_co2lab_state(model; pressure = p0))
    @test_throws ErrorException simulate_reservoir(case; info_level = -1)
end

# Sloping aquifer (sloping_aquifer.mat): geometry and well index against co2lab,
# and the well datum options on a short run.
function _ve_wells_vs_mrst(d)
    mesh = ve_mesh_from_mrst(d["G3"])
    cols, _ = ve_columns(mesh)
    ve = ve_domain(mesh, cols; perm_x = vec(d["rock"]["perm"]),
                   poro = vec(d["rock"]["poro"]))

    geo = ve.geometry

    rep = OracleReport("VE wells vs MRST co2lab")
    geomtol = Tolerance(1e-12, 1e-9)
    check!(rep, "caprock depth", geo.z_top, d["Gt_z"]; tol = geomtol)
    check!(rep, "H", geo.H, d["Gt_H"]; tol = geomtol)
    check!(rep, "dx (cellDims)", geo.dx, d["Gt_dx"]; tol = geomtol)
    check!(rep, "dy (cellDims)", geo.dy, d["Gt_dy"]; tol = geomtol)

    col, r = Int(d["well"]["cells"]), d["well"]["r"]
    last = Int(d["well_tpf"]["last_cell"])
    witol = Tolerance(1e-12, 1e-24)
    wi(c, ip) = ve_well_index(ve, c; radius = r, inner_product = ip)
    check!(rep, "WI, ip_simple", wi(col, :ip_simple), d["well"]["WI"]; tol = witol)
    check!(rep, "WI, ip_tpf", wi(col, :ip_tpf), d["well_tpf"]["WI"]; tol = witol)
    check!(rep, "WI, last column, ip_simple", wi(last, :ip_simple),
           d["well_tpf"]["last_WI_simple"]; tol = witol)
    check!(rep, "WI, last column, ip_tpf", wi(last, :ip_tpf),
           d["well_tpf"]["last_WI_tpf"]; tol = witol)

    passed(rep) || display(rep)
    @test passed(rep)

    _ve_well_datum_runs(d, ve, col, r)
end

# Rate-controlled injection with the datum at the caprock and at the bottom:
# the same flow, BHPs differing by the CO2 column rho_CO2·g·H.
function _ve_well_datum_runs(d, ve, col, r)
    geo, g, rhoc = ve.geometry, d["gravity"], d["rhoc"]
    fluid = co2lab_ve_fluid(krmax = (1.0, 1.0), co2_rho_ref = rhoc,
                            wat_rho_ref = d["rhow"], co2_mu_ref = d["muco2"],
                            wat_mu_ref = d["muw"])
    rc = (reference_pressure = d["p_ref"], compressibility = d["cf_rock"])

    p0 = ve_hydrostatic_pressure(geo; rho_w = d["rhow"], g = g)
    dt = vec(d["dT"])[1:8]

    function run(datum)
        w = setup_ve_well(ve, col; name = :inj, radius = r, datum = datum)
        model = setup_co2lab_model(ve.domain; wells = [w], fluid = fluid,
                                   rock_compressibility = rc)
        ctrl = co2_injector_control(TotalRateTarget(d["injection_rate"]), fluid)
        forces = setup_reservoir_forces(model, control = Dict(:inj => ctrl))

        prm = setup_parameters(model)
        # open-boundary stand-in (item 5): large end-column pore volumes
        prm[:Reservoir][:StaticFluidVolume][[1, end]] .*= 1e4

        case = JutulCase(model, dt, forces; parameters = prm,
                         state0 = setup_co2lab_state(model; pressure = p0))
        sim, cfg = setup_co2lab_simulator(case; info_level = -1, tol_cnv = 1e-9,
                                          tol_mb = 1e-12, max_nonlinear_iterations = 50)
        res = simulate_reservoir(case; simulator = sim, config = cfg)
        return prm[:inj][:PerforationGravityDifference][1], res
    end

    gdz_c, res_c = run(:caprock)
    gdz_b, res_b = run(:bottom)

    @test gdz_c == 0.0
    @test gdz_b ≈ -g * geo.H[col]

    sc, sb = res_c.states[end], res_b.states[end]
    @test maximum(sc[:sG]) > 0.01                  # CO2 has entered
    @test maximum(abs.(sc[:sG] .- sb[:sG])) < 1e-8
    @test maximum(abs.(sc[:Pressure] .- sb[:Pressure])) < 1e-2
    @test res_b.wells[:inj][:bhp][end] - res_c.wells[:inj][:bhp][end] ≈
          rhoc * g * geo.H[col] rtol = 1e-9

    # CO2 injected = CO2 in the aquifer + CO2 in the wellbore (its storage)
    injected = sum(res_c.wells[:inj][:mass_rate] .* dt)
    in_well = res_c.result.states[end][:inj][:TotalMasses][2, 1]
    @test sum(sc[:TotalMasses][2, :]) + in_well ≈ injected rtol = 1e-9
end

# ----------------------------------- Tests ------------------------------------

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

    @testset "VE columns: holes vs MRST" begin
        f = joinpath(refdir, "holes_case.mat")
        if !isfile(f)
            @info "Skipping: $f not generated"
        else
            d = matread(f)
            mesh = ve_mesh_from_mrst(d["G_in"])

            rep = OracleReport("VE columns, holes case")
            # :stack keeps MRST's discard=false cells, split at the 4 holes
            ck, dk = ve_columns(mesh; below_holes = :stack)
            check!(rep, ":stack: discarded", length(dk), 0; tol = exact)
            check!(rep, ":stack: columns.cells", ck.cells,
                   iv(d["keep"]["columns_cells"]); tol = exact)
            check!(rep, ":stack: columns", number_of_columns(ck),
                   length(iv(d["keep"]["columnPos"])) - 1 + 4; tol = exact)

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

            # opening sides: sealed fault sides only on request
            n_fault = count(sides.across_fault)
            @test length(ve_select_sides(sides, cols)) == length(sides) - n_fault
            @test length(ve_select_sides(sides, cols; include_across_fault = true)) ==
                  length(sides)
            right = ve_select_sides(sides, cols; directions = (:right,), j = 1:20)
            @test all(sides.direction[right] .== 2)
            @test all(cols.pillar[c][2] <= 20 for c in sides.column[right])
            @test_throws ArgumentError ve_select_sides(sides, cols; directions = (:top,))
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

            # initial state: brine hydrostatic at the caprock, as exampleVEJohansen.m:83
            s1 = matread(joinpath(refdir, "stage1_johansen.mat"))
            p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = 1020.0)
            p_mrst = 1020.0 * 9.80665 .* vec(s1["z_caprock"])
            @test maximum(abs.(p0 .- p_mrst) ./ p_mrst) < 1e-12

            # not the mid-column centroids, which lie H/2 lower on average
            zc = ve.domain[:cell_centroids][3, :]
            @test all(p0 .< 1020.0 * 9.80665 .* zc)

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

    @testset "stacked columns and connectivity (Stage 2 item 8)" begin
        layered(i, k) = 1000.0 + 10.0 * (k - 1)
        nc(m) = number_of_cells(m)
        pillar_ks(m, cols) = [[cell_ijk(m, c)[3] for c in column_cells(cols, col)]
                              for col in 1:number_of_columns(cols)]

        # a pinched-out layer in pillar 1 (zero thickness): its cells above and
        # below share a face, so no cut (co2lab would discard k = 3 there)
        pinched(i, k) = i == 1 && k >= 3 ? layered(i, k - 1) : layered(i, k)
        m = _cp_mesh(3, 3, pinched)
        cols, disc = ve_columns(m)
        @test isempty(disc) && pillar_ks(m, cols) == [[1, 3], [1, 2, 3], [1, 2, 3]]

        # a sealed face cuts: pillar 1 split between k = 1 and 2
        m = _cp_mesh(3, 3, layered)
        f = only(f for f in 1:number_of_faces(m) if
                 sort([cell_ijk(m, c) for c in m.faces.neighbors[f]]) ==
                 [(1, 1, 1), (1, 1, 2)])
        cols, disc = ve_columns(m; sealed_faces = [f])
        @test pillar_ks(m, cols) == [[1], [1, 2, 3], [1, 2, 3]] && length(disc) == 2
        cols, disc = ve_columns(m; sealed_faces = [f], below_holes = :stack)
        @test pillar_ks(m, cols)[1:2] == [[1], [2, 3]] && isempty(disc)

        # ...and the coarse face between the two stacked columns is sealed
        ve = ve_domain(m, cols; perm_x = fill(1e-13, nc(m)), poro = fill(0.2, nc(m)))
        N, T = ve.domain[:neighbors], ve.domain[:transmissibilities]
        @test [T[f] for f in axes(N, 2) if Set(N[:, f]) == Set([1, 2])] == [0.0]
        @test count(>(0), T) == size(N, 2) - 1

        # a shale missing from the grid in pillars 1-2 (k = 2), ending before 3
        full = _cp_mesh(3, 3, layered)
        shale = [c for c in 1:nc(full) if cell_ijk(full, c) in ((1, 1, 2), (2, 1, 2))]
        m = ve_trim(full, shale)

        cols, disc = ve_columns(m)
        @test pillar_ks(m, cols) == [[1], [1], [1, 2, 3]] && length(disc) == 2

        cols, disc = ve_columns(m; below_holes = :stack)
        @test pillar_ks(m, cols) == [[1], [3], [1], [3], [1, 2, 3]] && isempty(disc)
        @test cols.pillar == [(1, 1), (1, 1), (2, 1), (2, 1), (3, 1)]

        # column 5 (pillar 3) faces both stacked columns of pillar 2 on its left;
        # its left face at shale depth is exposed (no-flow), not an outer side
        cs = ve_column_sides(m, cols)
        left5 = findall(i -> cs.column[i] == 5 && cs.direction[i] == 1, 1:length(cs))
        @test cs.neighbour[left5] == [3, 4]
        bdir = JutulCO2Lab._boundary_directions(m)
        e = only(cs.exposed_faces)
        @test bdir[e] == 1 && cell_ijk(m, m.boundary_faces.neighbors[e]) == (3, 1, 2)
        ve = @test_logs ve_domain(m, cols; perm_x = fill(1e-13, nc(m)),
                                  poro = fill(0.2, nc(m)))
        @test ve.domain[:pillar_columns] == [2, 2, 2, 2, 1]
        @test all(ve.domain[:is_ve])
        @test Set(Set(ve.domain[:neighbors][:, f]) for f in 1:4) ==
              Set(Set.([[1, 3], [2, 4], [3, 5], [4, 5]]))

        sides, unassigned = ve_boundary_sides(m, cols, ve.domain)
        @test unassigned == cs.exposed_faces && !any(sides.across_fault)
        @test_throws ErrorException setup_co2lab_model(ve.domain)

        # layers that never communicate: warning
        cut = _cp_mesh(2, 3, layered; actnum = reshape([true, true, false, false,
                                                        true, true], 2, 1, 3))
        cols, _ = ve_columns(cut; below_holes = :stack)
        @test_logs (:warn, r"2 unconnected pieces") ve_domain(
            cut, cols; perm_x = fill(1e-13, 4), poro = fill(0.2, 4))

        # fine 3D cells are Stage 8; below_holes must be a known mode
        @test_throws ArgumentError ve_columns(m; is_ve = [true; fill(false, nc(m) - 1)])
        @test first(ve_columns(m; is_ve = fill(true, nc(m)))).pos == [1, 2, 3, 6]
        @test_throws ArgumentError ve_columns(m; below_holes = :bridge)

        # VE/fine branch points: compiled away for all-VE domains
        @test JutulCO2Lab._is_ve(nothing, 7)
        @test !JutulCO2Lab._is_ve([true, false], 2)
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

        # one-layer columns: VE column height = cell height
        domain[:column_height, Cells()] = fill(20.0, nx)

        inj = setup_well(domain, [nx - 3]; name = :inj)
        p0 = [1100.0 * 9.80665 * c[3] for c in eachcol(domain[:cell_centroids])]

        fluid = co2lab_ve_fluid(wat_rho_pvt = (4.3e-10, 100e5),
                                co2_rho_pvt = (4.3e-9, 100e5))
        rho_co2 = fluid.rho_ref[2]

        # same CO2 mass as 1e4 m³/day at 1.87 kg/m³; "surface" volumes are at the
        # reference density (760 kg/m³) in co2lab's convention
        q_in = 1e4 * 1.87 / 760.0 / day

        n_on, n_off = 10, 10
        dt = fill(10.0day, n_on + n_off)

        function run(; kw...)
            model = setup_co2lab_model(domain; wells = [inj], fluid = fluid, kw...)
            state0 = setup_co2lab_state(model; pressure = p0)

            on = InjectorControl(TotalRateTarget(q_in), [0.0, 1.0], density = rho_co2)
            f_on = setup_reservoir_forces(model, control = Dict(:inj => on))
            off = DisabledControl()
            f_off = setup_reservoir_forces(model, control = Dict(:inj => off))

            forces = vcat(fill(f_on, n_on), fill(f_off, n_off))

            case = JutulCase(model, dt, forces; state0 = state0,
                             parameters = setup_parameters(model))
            sim, cfg = setup_co2lab_simulator(case; info_level = -1)
            return model, simulate_reservoir(case; simulator = sim, config = cfg).states
        end

        injected = n_on * 10day * q_in * rho_co2
        co2_mass(s) = sum(s[:TotalMasses][2, :])

        @testset "configuration errors" begin
            # hysteresis: Stage 3; dissolution needs a rate; instantaneous: Stage 5
            hyst = co2lab_ve_fluid(residual = (0.2, 0.2))
            @test_throws ErrorException setup_co2lab_model(domain; fluid = hyst)
            @test_throws ErrorException setup_co2lab_model(domain; rs_sat = 0.03)
            @test_throws ErrorException setup_co2lab_model(domain; rs_sat = 0.03,
                                                           transfer_rate = Inf)
            @test_throws ErrorException co2lab_ve_fluid(co2_rho_pvt = (1, 2, 3, 4))

            # the VE fluid needs column heights
            plain = reservoir_domain(mesh, permeability = 0.5Darcy, porosity = 0.25)
            @test_throws ErrorException setup_co2lab_model(plain)

            # makeVEFluid options that belong to the grid or the model
            @test_throws ArgumentError co2lab_ve_fluid(rugosity = 1.0)
            @test_throws ArgumentError co2lab_ve_fluid(pvMult_fac = 1e-10)
        end

        @testset "rugosity (grid) and rock compressibility (model)" begin
            # rugosity per column from the domain; zero where there is none
            dom = reservoir_domain(mesh, permeability = 0.5Darcy, porosity = 0.25)
            dom[:column_height, Cells()] = fill(20.0, nx)
            dom[:rugosity, Cells()] = collect(1.0:nx) ./ 10
            prm(d; kw...) = setup_parameters(setup_co2lab_model(d; kw...))[:Reservoir]
            @test prm(dom)[:Rugosity] == collect(1.0:nx) ./ 10
            @test prm(domain)[:Rugosity] == zeros(nx)

            # compressible by default (makeVEFluid's values); static with nothing
            r = setup_co2lab_model(domain)[:Reservoir]
            @test r.secondary_variables[:FluidVolume] isa JutulCO2Lab.LinearPoreVolume
            @test r.secondary_variables[:FluidVolume].c == 1e-10
            r0 = setup_co2lab_model(domain; rock_compressibility = nothing)[:Reservoir]
            @test haskey(r0.parameters, :FluidVolume)
            @test !haskey(r0.secondary_variables, :FluidVolume)

            # ...and JutulDarcy's native form can be installed instead
            m = setup_co2lab_model(domain; rock_compressibility = nothing)
            JutulDarcy.set_rock_compressibility!(m; reference_pressure = 100e5,
                                                 compressibility = 1e-10)
            fv = m[:Reservoir].secondary_variables[:FluidVolume]
            @test fv isa JutulDarcy.LinearlyCompressiblePoreVolume

            # VE/fine branch points: no per-cell flag when every cell is VE
            S = r.secondary_variables
            @test S[:PlumeHeights].is_ve === nothing
            @test S[:RelativePermeabilities].is_ve === nothing
            @test S[:CapillaryPressure].is_ve === nothing
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

            # sGmax is updated after every converged sub-step, so it can exceed the
            # largest reported sG; it is never below it and never decreases
            runmax = accumulate((a, b) -> max.(a, b), [s[:sG] for s in st])
            @test all(all(s[:sGmax] .>= m .- 1e-12) for (s, m) in zip(st, runmax))
            @test all(all(b[:sGmax] .>= a[:sGmax]) for (a, b) in zip(st, st[2:end]))

            # the model's VE variables are the co2lab formulas (Stage 2 items 1, 3)
            fs = st[end]
            hh = JutulCO2Lab._plume_heights.(fs[:sG], fs[:sGmax], 20.0, 0.0, 0.0)
            @test fs[:PlumeHeights][1, :] ≈ first.(hh)
            @test fs[:PlumeHeights][2, :] ≈ last.(hh)

            # densities from the output pressure with co2lab's PVT (no dissolution)
            rho(ph) = fluid.rho_ref[ph] .*
                      JutulCO2Lab._linear.(1.0, Ref(fluid.rho_pvt[ph]), fs[:Pressure])
            @test fs[:CapillaryPressure][1, :] ≈
                  (rho(1) .- rho(2)) .* 9.80665 .* fs[:PlumeHeights][1, :]
            @test maximum(fs[:CapillaryPressure]) > 0

            # closed box: 0.5 % of the pore volume injected at c_t ≈ 5e-10/Pa raises
            # p by ≈ 100 bar; guard against runaway cases (1000s of bar)
            @test maximum(fs[:Pressure]) < 3 * maximum(p0)
        end

        @testset "dissolution on (rate-driven): unknowns p, sG, rs" begin
            rs_sat = 0.03
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

    @testset "VE wells (Stage 2 item 4)" begin
        wc = JutulCO2Lab._mrst_well_constant
        @test wc(10.0, 10.0, :ip_tpf) == 0.14
        @test wc(10.0, 10.0, :ip_simple) == 0.292

        # aspect ratio 2.5 rounds up, as MATLAB's round; 6 interpolates 5-8
        @test wc(25.0, 10.0, :ip_simple) == 0.262
        @test wc(10.0, 60.0, :ip_simple) ≈ 0.244 + (0.231 - 0.244) / 3
        @test_throws ArgumentError wc(1.0, 1.0, :ip_rt)

        geo1 = (z_top = [1000.0], H = [20.0])
        @test JutulCO2Lab._well_datum_depth(geo1, 1, :bottom) == 1020.0
        @test_throws ArgumentError JutulCO2Lab._well_datum_depth(geo1, 1, :top)

        fl = co2lab_ve_fluid(co2_rho_ref = 700)
        c = co2_injector_control(TotalRateTarget(1.0), fl)
        @test c.injection_mixture == [0.0, 1.0]
        @test c.mixture_density == 700.0
        @test c.phases == ((2, 1.0),)

        f = joinpath(refdir, "sloping_aquifer.mat")
        if !isfile(f)
            @info "Skipping: $f not generated (scripts/mrst_dump_sloping_aquifer.m)"
        else
            _ve_wells_vs_mrst(matread(f))
        end
    end

    @testset "VE boundary conditions (Stage 2 item 5)" begin
        f = joinpath(refdir, "ve_bc.mat")
        if !isfile(f)
            @info "Skipping: $f not generated (scripts/mrst_dump_ve_bc.m)"
        else
            _ve_bc_vs_mrst(matread(f))
        end

        f = joinpath(refdir, "sloping_aquifer.mat")
        if !isfile(f)
            @info "Skipping: $f not generated (scripts/mrst_dump_sloping_aquifer.m)"
        else
            _ve_bc_runs(matread(f))
        end
    end

    @testset "flow and well model vs MRST runs (Stage 2 item 7)" begin
        fs = [joinpath(refdir, f) for f in ("sloping_aquifer.mat", "mrst_ref_rate.mat",
                                            "mrst_ref_bhp.mat", "mrst_ref_src.mat")]
        if !all(isfile, fs)
            @info "Skipping: MRST runs not generated (scripts/stage2_mrst_comparison.jl)"
        else
            _ve_well_vs_mrst_runs(matread.(fs)...)
        end
    end

    @testset "VE reconstruction (Stage 2 item 6)" begin
        # one column of three cells, bottoms at 2, 5, 10 m below the caprock
        zb = [2.0, 5.0, 10.0]
        @test JutulCO2Lab._fill_degree(3.5, zb) == (1, 0.5)
        @test JutulCO2Lab._fill_degree(5.0, zb) == (2, 0.0)   # edge counts as full
        @test JutulCO2Lab._fill_degree(-1.0, zb) == (0, 0.0)
        @test JutulCO2Lab._fill_degree(12.0, zb) == (3, 0.0)

        # plume to 3.5 m, residual zone to 7.5 m: rw = 0.2, rg = 0.1
        s, smax = JutulCO2Lab._column_saturations(zb, 3.5, 7.5, 0.2, 0.1)
        @test s ≈ [0.8, 0.5 * 0.8 + 0.5 * 0.1, 0.5 * 0.1]
        @test smax ≈ [0.8, 0.8, 0.5 * 0.8]

        f = joinpath(refdir, "ve_reconstruction.mat")
        if !isfile(f)
            @info "Skipping: $f not generated (scripts/mrst_dump_ve_reconstruction.m)"
        else
            _ve_reconstruction_vs_mrst(matread(f))
        end
    end

    @testset "VE fluid vs MRST co2lab (Stage 2 items 1, 3)" begin
        f = joinpath(refdir, "ve_fluid.mat")
        if !isfile(f)
            @info "Skipping: $f not generated (scripts/mrst_dump_ve_fluid.m)"
        else
            d = matread(f)
            H = vec(d["H"]); g = d["gravity"]
            pvt = (co2_rho_ref = 760.0, wat_rho_ref = 1100.0,
                   co2_rho_pvt = (4.3e-9, 150e5), wat_rho_pvt = (4.3e-10, 100e5),
                   co2_mu_pvt = (2e-9, 150e5))
            # (fluid, caprock rugosity dh): rugosity is a grid property here
            cases = [(co2lab_ve_fluid(; residual = (0.0, 0.0), pvt...), 0.0),
                     (co2lab_ve_fluid(; residual = (0.27, 0.20), pvt...), 0.0),
                     (co2lab_ve_fluid(; residual = (0.27, 0.20), krmax = (0.9, 0.7),
                                      pvt...), 2.0)]
            # the dump's pore-volume multiplier: 1e-5 /bar about 100 bar
            pv_rock = JutulCO2Lab.LinearPoreVolume(100e5, 1e-10)
            tol = Tolerance(1e-12, 1e-14)
            rep = OracleReport("VE fluid vs MRST")
            for (c, (fl, dh)) in enumerate(cases)
                r = d["case$c"]
                # every sample (column) of every sample set, flattened
                sg, sgmax, p = vec(r["sg"]), vec(r["sgmax"]), vec(r["p"])
                Hs = repeat(H, size(r["sg"], 2))
                hh = JutulCO2Lab._plume_heights.(sg, sgmax, Hs, fl.res_water, fl.res_gas)
                h, hmax = first.(hh), last.(hh)
                kr = JutulCO2Lab._sharp_interface_relperm.(h, hmax, Hs, fl.krw, fl.krg,
                                                           dh)

                b(ph) = JutulCO2Lab._linear.(1.0, Ref(fl.rho_pvt[ph]), p)
                rho = [fl.rho_ref[ph] .* b(ph) for ph in 1:2]
                mu(ph) = JutulCO2Lab._linear.(fl.mu_ref[ph], Ref(fl.mu_pvt[ph]), p)

                check!(rep, "case $c: h", h, r["h"]; tol = tol)
                check!(rep, "case $c: h_max", hmax, r["hmax"]; tol = tol)
                check!(rep, "case $c: krW", first.(kr), r["krW"]; tol = tol)
                check!(rep, "case $c: krG", last.(kr), r["krG"]; tol = tol)
                check!(rep, "case $c: pcWG", (rho[1] .- rho[2]) .* g .* h, r["pcWG"];
                       tol = Tolerance(1e-12, 1e-9))
                check!(rep, "case $c: bW", b(1), r["bW"]; tol = tol)
                check!(rep, "case $c: bG", b(2), r["bG"]; tol = tol)
                check!(rep, "case $c: muW", mu(1), r["muW"]; tol = tol)
                check!(rep, "case $c: muG", mu(2), r["muG"]; tol = tol)

                pvm = 1 .+ pv_rock.c .* (p .- pv_rock.p_ref)
                check!(rep, "case $c: pore-volume multiplier", pvm, r["pvMult"];
                       tol = tol)
            end

            passed(rep) || display(rep)
            @test passed(rep)
        end

        # JutulDarcy adds pc to the non-reference phase: CO2 = brine + pcWG
        sys = co2lab_system()
        @test JutulDarcy.get_reference_phase_index(sys) == 1
        pc = [2.0e5 5.0e5]                         # one row, cells 1 and 2
        @test JutulDarcy.capillary_gradient(pc, 1, 2, 2, 1) == 3.0e5   # CO2
        @test JutulDarcy.capillary_gradient(pc, 1, 2, 1, 1) == 0.0     # brine
    end

end

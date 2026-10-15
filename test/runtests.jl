using Test
using JutulCO2Lab
using JutulCO2Lab.Oracle
using Jutul
using JutulDarcy
using MAT
using Statistics
using SparseArrays

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

# Transmissibility multipliers of `mesh`'s faces: `value` on `faces` (default
# 0: sealed), 1 elsewhere
function _multipliers(mesh, faces, value = 0.0)
    m = ones(number_of_faces(mesh))
    m[faces] .= value
    return m
end

# As _cp_mesh, but depths per pillar (p = 0..nx): layer k between pillars
# p and p+1 has corners top(p, k) and top(p + 1, k), so sloping layers are
# conforming (neighbouring cells share their faces fully)
function _cp_mesh_pillars(nx, nz, top)
    coord = Float64[]
    for j in 0:1, i in 0:nx
        append!(coord, [i * 100.0, j * 100.0, 0.0, i * 100.0, j * 100.0, 2000.0])
    end
    zc = [top(ii ÷ 2, k + s) for k in 1:nz for s in (0, 1)
          for jj in 1:2 for ii in 1:2nx]
    grid = Dict("cartDims" => [nx, 1, nz], "COORD" => coord, "ZCORN" => zc,
                "ACTNUM" => fill(true, nx, 1, nz))
    return JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
end

# --------------------- Helpers: a stand-in vertical model ---------------------

# The sharp interface's formulas as a model of its own, so that it runs through
# the generic path that later vertical models use (densities as an input)
struct _StandInSharp <: JutulCO2Lab.VerticalModel end
JutulCO2Lab._column_heights(::_StandInSharp, i, sG, sGmax, H, drho, rw, rg) =
    JutulCO2Lab._plume_heights(sG, sGmax, H, rw, rg)
JutulCO2Lab._column_relperm(::_StandInSharp, i, sG, sGmax, h, h_max, H, drho, dh, krw,
                            krg, rw, rg) =
    JutulCO2Lab._sharp_interface_relperm(h, h_max, H, krw, krg, dh)

# --------------- Helpers: comparisons with MRST dumps and runs ----------------

include("atlas_checks.jl")
include("trap_checks.jl")
include("inventory_checks.jl")
include("capacity_checks.jl")
include("fringe_checks.jl")

# The VE sloping-aquifer domain of an MRST dump (3D grid G3 and rock).
function _ve_case(d)
    mesh = ve_mesh_from_mrst(d["G3"])
    cols, _ = ve_columns(mesh)
    ve = ve_domain(mesh, cols; perm_x = vec(d["rock"]["perm"]),
                   poro = vec(d["rock"]["poro"]))

    sides, _ = ve_boundary_sides(mesh, cols, ve.domain)
    return ve, cols, sides
end

# capfringe_si.mat: the sharp-interface half
# of co2lab's capfringeFigures.m — relperm and pcWG over sG without and with
# hysteresis, and the reconstructed profile of one state, on its single column.
function _capfringe_si_vs_mrst(d)
    mesh = ve_mesh_from_mrst(d["G3"])
    cols, _ = ve_columns(mesh)
    ve = ve_domain(mesh, cols; perm_x = fill(1e-13, 1000), poro = fill(0.3, 1000))
    H, rw, rg = ve.geometry.H[1], d["srw"], d["src"]
    fl = co2lab_ve_fluid(default_pvt = nothing,
                         residual = (rw, rg), krmax = Tuple(vec(d["krmax"])),
                         wat_rho_ref = d["rhoW"], co2_rho_ref = d["rhoG"])

    # the model's plume heights, relperm and pcWG for (sG, sGmax)
    function curves(s, smax)
        hh = JutulCO2Lab._plume_heights.(s, max.(min.(smax, 1 - rw), s), H, rw, rg)
        kr = JutulCO2Lab._sharp_interface_relperm.(first.(hh), last.(hh), H, fl.krw,
                                                   fl.krg, 0.0)
        pc = (d["rhoW"] - d["rhoG"]) * d["gravity"] .* first.(hh)
        return first.(kr), last.(kr), pc
    end

    rep = OracleReport("capfringeFigures.m, sharp interface, vs MRST")
    tol, tol_pc = Tolerance(1e-12, 1e-14), Tolerance(1e-12, 1e-9)
    s = vec(d["s"])
    cases = (("", s, ""), (" with hysteresis", vec(d["smax_hyst"]), "_hyst"))
    for (label, smax, sfx) in cases
        krw, krg, pc = curves(s, smax)
        check!(rep, "krW" * label, krw, vec(d["krW" * sfx]); tol = tol)
        check!(rep, "krG" * label, krg, vec(d["krG" * sfx]); tol = tol)
        check!(rep, "pcWG" * label, pc, vec(d["pcWG" * sfx]); tol = tol_pc)
    end

    st = Dict(:Pressure => [d["pressure"]], :sG => [d["S"]], :sGmax => [d["Smax"]])
    o = ve_reconstruct(st, ve, cols, fl)
    check!(rep, "profile: sG", o.sG, vec(d["s3D"]); tol = Tolerance(0.0, 1e-12))
    check!(rep, "profile: sGmax", o.sGmax, vec(d["smax3D"]); tol = Tolerance(0.0, 1e-12))
    return rep
end

# mrst_ref_{rate,bhp,src}.mat: co2lab's sloping aquifer run by MRST on Julia's
# geometry, same report steps, tight tolerances. With a source term instead of
# the well, Julia matches to solver precision (measured 5e-10); with the well,
# to the wellbore-storage difference (measured BHP 1e-5, rate 3e-4, pressure
# 5e-6, sG 1.3e-5).
function _ve_well_vs_mrst_runs(d, ref_rate, ref_bhp, ref_src)
    ve, cols, sides = _ve_case(d)
    g, rho_w = d["gravity"], d["rhow"]
    fluid = co2lab_ve_fluid(default_pvt = nothing,
                            krmax = (1.0, 1.0), co2_rho_ref = d["rhoc"],
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
                         state0 = setup_co2lab_state(model; pressure = p0,
                                                     forces = f_on),
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
        fl = co2lab_ve_fluid(default_pvt = nothing,
                             residual = res, co2_rho_ref = 700, wat_rho_ref = 1050,
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

    # our half transmissibility, from the 3D grid
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
    fluid = co2lab_ve_fluid(default_pvt = nothing,
                            krmax = (1.0, 1.0), co2_rho_ref = rhoc, wat_rho_ref = rhow,
                            co2_mu_ref = d["muco2"], wat_mu_ref = d["muw"])
    rc = (reference_pressure = d["p_ref"], compressibility = d["cf_rock"])

    p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = rhow, g = g)
    hydro(z) = rhow * g * z
    all_open = ve_boundary_conditions(ve, sides, ve_select_sides(sides, cols);
                                      pressure = hydro)

    dt = vec(d["dT"])[1:8]

    function run(bcs, wells, control; fl = fluid)
        model = setup_co2lab_model(ve.domain; wells = wells, fluid = fl,
                                   rock_compressibility = rc)
        push!(model[:Reservoir].output_variables, :RelativePermeabilities,
              :PhaseViscosities, :PhaseMassDensities, :ShrinkageFactors)
        forces = setup_reservoir_forces(model; control = control, bc = bcs)
        case = JutulCase(model, dt, forces; parameters = setup_parameters(model),
                         state0 = setup_co2lab_state(model; pressure = p0,
                                                     forces = forces))
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

    # the wellbore starts filled with CO2 and stays so: count the reservoir only
    last = res.states[end]
    total = sum(last[:Reservoir][:TotalMasses][2, :])

    # brine is pushed out; the CO2 has not reached the boundary yet
    @test exited[1] > 0
    @test (total + exited[2]) / injected ≈ 1 rtol = 1e-5

    # with dissolution: CO2 is conserved when brine carries
    # its dissolved CO2 across (the default); co2lab's booking has no
    # dissolved part in the exiting CO2
    fl_dis = co2lab_ve_fluid(default_pvt = nothing,
                             krmax = (1.0, 1.0), co2_rho_ref = rhoc, wat_rho_ref = rhow,
                             co2_mu_ref = d["muco2"], wat_mu_ref = d["muw"],
                             dissolution = true, dis_rate = 3e-11, dis_max = 0.05)
    for mode in (:carry, :co2lab)
        bcs_d = ve_boundary_conditions(ve, sides, right; pressure = hydro,
                                       dissolved = mode)
        res = run(bcs_d, [w], Dict(:inj => ctrl); fl = fl_dis)
        sub, dts, _ = Jutul.expand_to_ministeps(res)
        q = [ve_boundary_mass_fluxes(bcs_d, s[:Reservoir]) for s in sub]
        injected = sum(rate(s) * t for (s, t) in zip(sub, dts))
        exited_co2 = sum(sum(qk[2, :]) * t for (qk, t) in zip(q, dts))
        last = res.states[end]
        total = sum(last[:Reservoir][:TotalMasses][2, :])
        @test (total + exited_co2) / injected ≈ 1 rtol = 1e-5
        mode == :co2lab && @test all(qk -> all(iszero, qk[3, :]), q)
    end
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
    fluid = co2lab_ve_fluid(default_pvt = nothing,
                            krmax = (1.0, 1.0), co2_rho_ref = rhoc,
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

        @testset "grid tolerances (STAGE1_TOL) are defined for every quantity" begin
            for k in (:height, :pore_volume, :transmissibility, :trans_mult, :well_index)
                @test haskey(STAGE1_TOL, k)
            end
        end
    end

    # Comparisons with MRST reference dumps (.mat), from the directory in
    # JUTULCO2LAB_REFERENCE_DATA (default: analysis/reference-data next to the
    # package's directory); each skips when its dump is missing.
    refdir = get(ENV, "JUTULCO2LAB_REFERENCE_DATA",
                 joinpath(@__DIR__, "..", "..", "analysis", "reference-data"))
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
            @test all(cols.column_ij[c][2] <= 20 for c in sides.column[right])
            @test_throws ArgumentError ve_select_sides(sides, cols; directions = (:top,))
        end
    end

    @testset "sloping aquifer built in Julia vs MRST" begin
        # co2lab's ramp-up: 7 halvings, then whole steps (MRST's rampupTimesteps)
        dts = rampup_timesteps(24.0, 1.0, 7)
        @test dts == [1 / 128; 1 ./ 2.0 .^ (7:-1:1); fill(1.0, 23)]
        f = joinpath(refdir, "sloping_aquifer.mat")
        if !isfile(f)
            @info "Skipping: $f not generated"
        else
            d = matread(f)
            aq = sloping_aquifer(; nz = 5, co2_rho = d["rhoc"], co2_mu = d["muco2"])
            m0 = ve_mesh_from_mrst(d["G3"])
            # make_testgrid's nodes exactly, the cells in MRST's order
            X0 = Set(Tuple(r) for r in eachrow(d["G3"]["nodes"]["coords"]))
            @test Set(Tuple(x) for x in aq.mesh.node_points) == X0
            @test [cell_ijk(aq.mesh, c) for c in 1:number_of_cells(aq.mesh)] ==
                  [cell_ijk(m0, c) for c in 1:number_of_cells(m0)]
            @test isapprox(aq.ve.geometry.H, vec(d["Gt_H"]); rtol = 1e-12)
            @test isapprox(aq.p_ref, d["p_ref"]; rtol = 1e-15)
            @test aq.well_column == Int(d["well"]["cells"])
            @test aq.mass_rate / aq.co2_rho == d["injection_rate"]
            @test aq.dt == vec(d["dT"]) && aq.injecting == (vec(d["step_control"]) .== 1)
            @test aq.sides.column[aq.open_sides] == [Int(d["bc"]["cell"])]
        end
    end

    @testset "hydrostatic brine of the fluid (ve_hydrostatic_profile)" begin
        aq = sloping_aquifer(; nz = 5, co2_rho = 700.0, co2_mu = 5e-5)
        geo = aq.ve.geometry
        # constant brine: co2lab's rho_w·g·z
        f0 = co2lab_ve_fluid(; default_pvt = nothing, wat_rho_ref = 1050.0)
        p = ve_hydrostatic_profile(f0, aq.ve).(geo.z_top)
        @test isapprox(p, ve_hydrostatic_pressure(geo; rho_w = 1050.0); rtol = 1e-12)
        # JutulDarcy's brine (30 °C here): at rest between columns, which
        # 1050 kg/m³ is not
        f1 = co2lab_ve_fluid()
        N, gdz = aq.ve.domain[:neighbors], aq.ve.domain[:gdz]
        function imbalance(p)
            rho = [JutulCO2Lab._phase_value(f1.rho_ref[1], f1.rho_pvt[1], q, 303.15)
                   for q in p]
            return maximum(abs((p[N[2, k]] - p[N[1, k]]) +
                               (rho[N[1, k]] + rho[N[2, k]]) / 2 * gdz[k]) /
                           (1000 * abs(gdz[k]) + 1e-12) for k in axes(N, 2))
        end
        @test imbalance(ve_hydrostatic_profile(f1, aq.ve).(geo.z_top)) < 1e-5
        @test imbalance(ve_hydrostatic_pressure(geo; rho_w = 1050.0)) > 1e-2
    end

    @testset "Johansen VE domain from the data set vs MRST" begin
        # without the data: an error naming the download address
        withenv("CO2LAB_ATLAS_DATA" => mktempdir()) do
            err = try johansen_grid_section(); nothing catch e; e end
            @test err isa ArgumentError && occursin("johansen.zip", err.msg)
        end
        fg = joinpath(refdir, "johansen_grids.mat")
        dir = johansen_data_directory()
        if !(isfile(fg) && isfile(joinpath(dir, "NPD5_Permeability.txt")))
            @info "Skipping: Johansen data set or reference data not present"
        else
            jo = johansen_ve_domain()
            d = matread(fg)
            @test jo.grid["layers"] == 6:10
            # the same cells in the same order, the same rock bit for bit
            G3 = d["G3"]
            nx, ny = Int.(G3["cartDims"][1:2])
            ijk0 = [((m - 1) % nx + 1, ((m - 1) ÷ nx) % ny + 1, (m - 1) ÷ (nx * ny) + 1)
                    for m in Int.(vec(G3["cells"]["indexMap"]))]
            ijk = [cell_ijk(jo.mesh, c) .+ (0, 0, 5) for c in 1:number_of_cells(jo.mesh)]
            @test ijk == ijk0
            at(key) = [jo.grid[key][c[1], c[2], c[3] - 5] for c in ijk0]
            K = d["rock3D"]["perm"]
            @test at("PERMX") == K[:, 1] && at("PERMY") == K[:, 2] &&
                  at("PERMZ") == K[:, 3] && at("PORO") == vec(d["rock3D"]["poro"])
            # the same columns; one side per outer edge of co2lab's Gt (GeoEnergyIO's
            # side directions corrected); co2lab's open boundary from its ranges
            mcols, _ = ve_columns(ve_mesh_from_mrst(G3))
            @test jo.columns.column_ij == mcols.column_ij
            Nt = Int.(d["Gt"]["faces"]["neighbors"])
            @test length(jo.sides) == count((Nt[:, 1] .== 0) .| (Nt[:, 2] .== 0))
            sel(dir, i, j) = ve_select_sides(jo.sides, jo.columns; directions = (dir,),
                                             i, j)
            n_open = length(union(sel(:lower, 1:94, 1:4), sel(:left, 1:20, nothing),
                                  sel(:right, nothing, 90:100),
                                  sel(:upper, 1:42, 50:100)))
            @test n_open == 234                         # co2lab's bcIxVE
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
            fault_i = Set(cols.column_ij[c][1] for c in sides.column[sides.across_fault])
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

            # The fault's two exposed strips, no-flow: above i = 3, the upper part
            # of the i = 2 columns' right side (cells k = 1, 2, tagged :right); below
            # i = 2, the lower part of the i = 3 columns' left side (cells k = 2, 3,
            # tagged :left; GeoEnergyIO 1.1.34 tagged them :right, fixed in the fork,
            # sintefmath/GeoEnergyIO.jl #66)
            bdir = JutulCO2Lab._boundary_directions(mesh)
            i_of(f) = cell_ijk(mesh, mesh.boundary_faces.neighbors[f])[1]
            @test length(unassigned) == 8
            @test sort([(i_of(f), bdir[f]) for f in unassigned]) ==
                  sort([fill((2, 2), 4); fill((3, 1), 4)])
        end

        dir = JutulCO2Lab._logical_direction
        @test_throws ErrorException dir((1, 1, 1), (2, 2, 1))   # diagonal
        @test dir((3, 4, 5), (3, 4, 2)) == 5                    # top
        @test dir((3, 4, 5), (4, 4, 9)) == 2                    # i+, across a fault
    end

    @testset "stacked columns and connectivity" begin
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
        cols, disc = @test_logs (:warn, r"below_holes = :discard") ve_columns(
            m; transmissibility_multiplier = _multipliers(m, [f]))
        @test pillar_ks(m, cols) == [[1], [1, 2, 3], [1, 2, 3]] && length(disc) == 2
        cols, disc = ve_columns(m; transmissibility_multiplier = _multipliers(m, [f]),
                                below_holes = :stack)
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
        @test cols.column_ij == [(1, 1), (1, 1), (2, 1), (2, 1), (3, 1)]

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
        # stacked blocks: their faces take the transition flux
        law = setup_co2lab_model(ve.domain)[:Reservoir].equations[:mass_conservation]
        @test Jutul.flux_type(law) isa JutulCO2Lab.VEFlux

        # a well over pillar 1's two stacked columns: one perforation at each
        # caprock, each with its column's index; the wellbore's hydrostatic
        # term from the datum (the upper caprock, or the lower column's bottom)
        g = Jutul.gravity_constant
        gdz(datum) = setup_parameters(setup_co2lab_model(ve.domain; wells = [
            setup_ve_well(ve, [1, 2]; name = :inj, datum)]))[:inj]
        @test gdz(:caprock)[:WellIndices] ≈ ve_well_index.(Ref(ve), [1, 2])
        @test gdz(:caprock)[:PerforationGravityDifference] ≈ [0.0, 20g]
        @test gdz(:bottom)[:PerforationGravityDifference] ≈ [-30g, -10g]

        # layers that never communicate: warning
        cut = _cp_mesh(2, 3, layered; actnum = reshape([true, true, false, false,
                                                        true, true], 2, 1, 3))
        cols, _ = ve_columns(cut; below_holes = :stack)
        @test_logs (:warn, r"2 unconnected pieces") ve_domain(
            cut, cols; perm_x = fill(1e-13, 4), poro = fill(0.2, 4))

        # fine 3D cells are blocks of their own; below_holes must be
        # a known mode
        hyb, _ = ve_columns(m; is_ve = [true; fill(false, nc(m) - 1)])
        @test count(hyb.is_ve) == 1 &&
              all(length(column_cells(hyb, b)) == 1 for b in findall(!, hyb.is_ve))
        @test first(ve_columns(m; is_ve = fill(true, nc(m)))).pos == [1, 2, 3, 6]
        @test_throws ArgumentError ve_columns(m; below_holes = :bridge)

        # VE/fine branch points: compiled away for all-VE domains
        @test JutulCO2Lab._is_ve(nothing, 7)
        @test !JutulCO2Lab._is_ve([true, false], 2)
    end

    @testset "net-to-gross: net rock, cells without net rock removed" begin
        layered(i, k) = 1000.0 + 10.0 * (k - 1)
        full = _cp_mesh(3, 3, layered)
        n = number_of_cells(full)
        ijk = [cell_ijk(full, c) for c in 1:n]
        # net-to-gross 0 in a shale (pillars 1-2, k = 2), 0.2-1 elsewhere
        ntg = [i <= 2 && k == 2 ? 0.0 : 0.2 * (i + k - 1) for (i, _, k) in ijk]
        poro = [0.1 + 0.05k for (_, _, k) in ijk]
        perm = fill(1e-13, n)
        nk(m, cols) = [[cell_ijk(m, c)[3] for c in column_cells(cols, b)]
                       for b in 1:number_of_columns(cols)]

        # the shale removed: a hole, so stacked columns above and below it
        net = net_rock(full; poro, perm_x = perm, ntg)
        keep = findall(>(0), ntg)
        @test net.cells == keep && number_of_cells(net.mesh) == length(keep)
        @test net.poro ≈ poro[keep] .* ntg[keep] && net.perm_x ≈ perm[keep] .* ntg[keep]
        @test isnothing(net.perm_y) && isnothing(net.perm_z)
        cols, _ = ve_columns(net.mesh; below_holes = :stack)
        @test nk(net.mesh, cols) == [[1], [3], [1], [3], [1, 2, 3]]

        # a VE column's pore volume: the sum of the cells' vol·poro·ntg (not
        # the product of column averages)
        ve = ve_domain(net.mesh, cols; perm_x = net.perm_x, poro = net.poro)
        vol = reservoir_domain(full; permeability = perm, porosity = poro)[:volumes]
        c3 = findall(c -> c[1] == 3, ijk)
        @test JutulDarcy.pore_volume(ve.domain)[5] ≈ sum((vol .* poro .* ntg)[c3])
        @test !haskey(ve.domain, :net_to_gross)

        # with perm_z given: JutulDarcy's net_to_gross (ECLIPSE's NTG), the
        # same pore volumes and transmissibilities on every face
        pos, kz = max.(ntg, 0.1), fill(2e-14, n)
        net = net_rock(full; poro, perm_x = perm, perm_z = kz, ntg = pos)
        @test net.mesh === full && net.perm_z == kz
        ours = reservoir_domain(full; porosity = net.poro,
                                permeability = [net.perm_x'; net.perm_x'; net.perm_z'])
        jd = reservoir_domain(full; porosity = poro, permeability = [perm'; perm'; kz'],
                              net_to_gross = pos)
        T(d) = JutulDarcy.reservoir_transmissibility(d)
        @test JutulDarcy.pore_volume(ours) ≈ JutulDarcy.pore_volume(jd)
        @test T(ours) ≈ T(jd)

        # one permeability per cell: scaled in all directions, so vertical
        # faces (within a pillar) carry less than with JutulDarcy's
        one = net_rock(full; poro, perm_x = perm, ntg = pos)
        T1 = T(reservoir_domain(full; porosity = one.poro, permeability = one.perm_x))
        Tj = T(reservoir_domain(full; porosity = poro, permeability = perm,
                                net_to_gross = pos))
        N = get_neighborship(full)
        vert = [ijk[N[1, f]][1] == ijk[N[2, f]][1] for f in axes(N, 2)]
        @test T1[.!vert] ≈ Tj[.!vert] && all(T1[vert] .< Tj[vert])

        # a threshold; no net-to-gross; bad values
        @test net_rock(full; poro, perm_x = perm, ntg, ntg_min = 0.5).cells ==
              findall(>(0.5), ntg)
        none = net_rock(full; poro, perm_x = perm, ntg = nothing)
        @test none.mesh === full && none.poro == poro && none.poro !== poro
        @test_throws ArgumentError net_rock(full; poro, perm_x = perm, ntg = ntg[2:end])
        @test_throws ArgumentError net_rock(full; poro, perm_x = perm,
                                            ntg = [NaN; ntg[2:end]])
        @test_throws ArgumentError net_rock(full; poro, perm_x = perm, ntg = ntg .- 0.1)
    end

    @testset "barriers: transmissibility multipliers, low vertical permeability" begin
        J = JutulCO2Lab
        layered(i, k) = 1000.0 + 10.0 * (k - 1)
        mesh = _cp_mesh(3, 3, layered)
        nc, nf = number_of_cells(mesh), number_of_faces(mesh)
        ijk = [cell_ijk(mesh, c) for c in 1:nc]
        N = mesh.faces.neighbors
        nk(cols) = [[ijk[c][3] for c in column_cells(cols, b)]
                    for b in 1:number_of_columns(cols)]
        face(a, b) = only(f for f in 1:nf if Set(ijk[[N[f]...]]) == Set([a, b]))
        perm, poro = fill(1e-13, nc), fill(0.2, nc)
        Tf = JutulDarcy.reservoir_transmissibility(
            reservoir_domain(mesh; permeability = perm, porosity = poro))

        # between k = 1 and 2: pillar 1 sealed (0), pillar 2 a barrier (0.01
        # below the threshold 0.1), pillar 3 not (0.5)
        m = ones(nf)
        f1, f2, f3 = (face((i, 1, 1), (i, 1, 2)) for i in 1:3)
        m[[f1, f2, f3]] = [0.0, 0.01, 0.5]
        cols, _ = ve_columns(mesh; transmissibility_multiplier = m,
                             barrier_multiplier = 0.1, below_holes = :stack)
        @test nk(cols) == [[1], [2, 3], [1], [2, 3], [1, 2, 3]]
        @test cols.sealed == [f1] && cols.leaky == [f2] && cols.multiplier == m

        # the barrier's transmissibility: the fine face's times its multiplier
        d = ve_domain(mesh, cols; perm_x = perm, poro).domain
        T, cat, Nc = d[:transmissibilities], d[:face_category], d[:neighbors]
        across(a, b) = only(f for f in axes(Nc, 2) if Set(Nc[:, f]) == Set([a, b]))
        @test cat[across(3, 4)] == J._FACE_VE_VERTICAL && T[across(3, 4)] ≈ 0.01Tf[f2]
        @test cat[across(1, 2)] == J._FACE_SEALED && T[across(1, 2)] == 0

        # without a threshold only the multiplier 0 cuts
        cols, _ = ve_columns(mesh; transmissibility_multiplier = m, below_holes = :stack)
        @test nk(cols) == [[1], [2, 3], [1, 2, 3], [1, 2, 3]]

        # pure VE, a lateral multiplier: the continuous VE–VE face's coarse
        # transmissibility times the multiplier
        lat = [f for f in 1:nf if Set(ijk[n][1] for n in N[f]) == Set([1, 2])]
        ml = ones(nf)
        ml[lat] .= 0.5
        T0 = JutulDarcy.reservoir_transmissibility(
            ve_domain(mesh, first(ve_columns(mesh)); perm_x = perm, poro).domain)
        dl = ve_domain(mesh, first(ve_columns(mesh; transmissibility_multiplier = ml));
                       perm_x = perm, poro).domain
        @test all(==(J._FACE_VE), dl[:face_category])
        half = [Set(dl[:neighbors][:, f]) == Set([1, 2]) for f in eachindex(T0)]
        @test dl[:transmissibilities] ≈ T0 .* ifelse.(half, 0.5, 1.0) && count(half) == 1

        # a cell of low vertical permeability: a fine 3D cell
        kz = copy(perm)
        kz[findfirst(==((1, 1, 2)), ijk)] = 1e-18
        cols, _ = ve_columns(mesh; perm_z = kz, barrier_perm_z = 1e-16)
        @test nk(cols)[1:3] == [[1], [2], [3]] && cols.is_ve[1:3] == [true, false, true]

        # thresholds need their data; :discard with barriers warns
        @test_throws ArgumentError ve_columns(mesh; barrier_multiplier = 0.1)
        @test_throws ArgumentError ve_columns(mesh; barrier_perm_z = 1e-16)
        @test_throws ArgumentError ve_columns(mesh; transmissibility_multiplier = m[2:end])
        @test_logs (:warn, r"below_holes = :discard") ve_columns(
            mesh; transmissibility_multiplier = m, barrier_multiplier = 0.1)
    end

    @testset "hybrid grid: blocks and face categories" begin
        J = JutulCO2Lab
        nc(m) = number_of_cells(m)
        dom(m, cols; ff = :virtual) = ve_domain(m, cols; perm_x = fill(1e-13, nc(m)),
                                                poro = fill(0.2, nc(m)), fault_flux = ff)
        counts(d) = Dict(c => count(==(c), d[:face_category])
                         for c in unique(d[:face_category]))
        # a conforming sloping layer: every lateral face continuous, nothing
        # changes for pure VE (no transmissibilities set)
        slope = _cp_mesh_pillars(4, 4, (p, k) -> 1000.0 + 3.0 * p + 2.0 * (k - 1))
        cols, _ = ve_columns(slope)
        d = dom(slope, cols).domain
        @test counts(d) == Dict(J._FACE_VE => 3) && !haskey(d, :transmissibilities)

        # a fault offset by 3 m (cells 2 m): the top cells do not connect —
        # discontinuous with :virtual, continuous with co2lab's flux
        fault = _cp_mesh(4, 4, (i, k) -> 1000.0 + (i > 2 ? 3.0 : 0.0) + 2.0 * (k - 1))
        cols, _ = ve_columns(fault)
        @test counts(dom(fault, cols; ff = :co2lab).domain) == Dict(J._FACE_VE => 3)
        d = dom(fault, cols).domain
        f = only(findall(==(J._FACE_VE_DISCONTINUOUS), d[:face_category]))
        @test Set(d[:neighbors][:, f]) == Set([2, 3])
        # each side's cells at the face, below its own top: 2–8 m and 0–6 m
        @test d[:face_intervals][:, f] ≈ (d[:neighbors][1, f] == 2 ? [2, 8, 0, 6] :
                                          [0, 6, 2, 8])
        # the transmissibility summed from the fine grid's
        fine = reservoir_domain(fault; permeability = 1e-13, porosity = 0.2)
        Tf = JutulDarcy.reservoir_transmissibility(fine)
        cf = ve_domain(fault, cols; perm_x = fill(1e-13, nc(fault)),
                       poro = fill(0.2, nc(fault)), fault_flux = :virtual).coarse
        @test d[:transmissibilities][f] ≈ sum(Tf[cf.coarse_faces_to_fine[f]])

        # whole fine pillars, given as (i, j) pairs or as a function of them
        ijk = [cell_ijk(slope, c) for c in 1:nc(slope)]
        @test ve_except_pillars(slope, [(1, 1), (3, 1)]) ==
              ve_except_pillars(slope, (i, j) -> i in (1, 3)) ==
              [!(x[1] in (1, 3)) for x in ijk]
        cols, _ = ve_columns(slope; is_ve = ve_except_pillars(slope, [(2, 1)]))
        @test cols.is_ve == [true, false, false, false, false, true, true]

        # fine cells: pillar 1 and the top two cells of pillar 3
        isve = [!(x[1] == 1 || (x[1] == 3 && x[3] <= 2)) for x in ijk]
        cols, _ = ve_columns(slope; is_ve = isve)
        @test cols.is_ve == [false, false, false, false, true, false, false, true, true]
        ve = dom(slope, cols)
        d = ve.domain
        # fine–fine 4, VE–fine 9 (lateral and the one vertical), and the VE
        # blocks beside pillar 3's lower block discontinuous (2)
        @test counts(d) == Dict(J._FACE_FINE => 4, J._FACE_VE_FINE => 9,
                                J._FACE_VE_DISCONTINUOUS => 2)
        @test d[:is_ve] == cols.is_ve
        # the vertical VE–fine face: pillar 3's fine k = 2 above its VE block
        fv = only(f for f in axes(d[:neighbors], 2) if
                  Set(d[:neighbors][:, f]) == Set([7, 8]))
        side = d[:neighbors][1, fv] == 8 ? 1 : 2
        # the VE block's top cell (0–2 m below its top); common depth the fine
        # cell's centroid
        @test d[:face_intervals][2side - 1:2side, fv] ≈ [0.0, 2.0]
        @test d[:face_depth][fv] ≈ d[:cell_centroids][3, 7]
        # fine cells' gravity from their centroids
        ff = only(f for f in axes(d[:neighbors], 2) if
                  Set(d[:neighbors][:, f]) == Set([1, 2]))
        @test abs(d[:gdz][ff]) ≈ Jutul.gravity_constant *
                                 abs(d[:cell_centroids][3, 2] - d[:cell_centroids][3, 1])
        @test_throws ErrorException setup_co2lab_model(d)

        # a layer between k = 2 and 3: a barrier (multiplier 0.01) — stacked
        # blocks, vertical transitions with the fine grid's transmissibility
        # times 0.01; sealed (0) — no flow
        N = slope.faces.neighbors
        lay = [f for f in 1:number_of_faces(slope) if
               ijk[N[f][1]][1:2] == ijk[N[f][2]][1:2] &&
               minmax(ijk[N[f][1]][3], ijk[N[f][2]][3]) == (2, 3)]
        cols, _ = ve_columns(slope; transmissibility_multiplier = _multipliers(slope,
                             lay, 0.01), barrier_multiplier = 0.1, below_holes = :stack)
        d = dom(slope, cols).domain
        @test count(==(J._FACE_VE_VERTICAL), d[:face_category]) == 4
        @test all(d[:transmissibilities][d[:face_category] .== J._FACE_VE_VERTICAL] .> 0)
        cols, _ = ve_columns(slope; transmissibility_multiplier = _multipliers(slope,
                             lay), below_holes = :stack)
        d = dom(slope, cols).domain
        @test count(==(J._FACE_SEALED), d[:face_category]) == 4
        @test all(d[:transmissibilities][d[:face_category] .== J._FACE_SEALED] .== 0)
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

    @testset "co2lab fluid system" begin
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

        fluid = co2lab_ve_fluid(default_pvt = nothing, wat_rho_pvt = (4.3e-10, 100e5),
                                co2_rho_pvt = (4.3e-9, 100e5))
        rho_co2 = fluid.rho_ref[2]

        # same CO2 mass as 1e4 m³/day at 1.87 kg/m³; "surface" volumes are at the
        # reference density (760 kg/m³) in co2lab's convention
        q_in = 1e4 * 1.87 / 760.0 / day

        n_on, n_off = 10, 10
        dt = fill(10.0day, n_on + n_off)

        function run(; fluid = fluid, kw...)
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
            # unknowns by configuration: hysteresis with rate-driven
            # dissolution couples the two equations; instantaneous dissolution
            # (dis_rate 0 or Inf) has co2lab's switching unknown instead of sG, Rs
            hd = co2lab_ve_fluid(default_pvt = nothing,
                                 residual = (0.2, 0.2), dissolution = true,
                                 dis_rate = 1e-8)
            r = setup_co2lab_model(domain; fluid = hd)[:Reservoir]
            @test collect(keys(r.primary_variables)) == [:Pressure, :sG, :Rs, :sGmax]
            @test r.equations[:hysteresis].dissolution === r.equations[:dissolved_co2]
            for rate in (0.0, Inf)
                inst = co2lab_ve_fluid(default_pvt = nothing,
                                       residual = (0.2, 0.2), dissolution = true,
                                       dis_rate = rate)
                @test instantaneous_dissolution(inst)
                r = setup_co2lab_model(domain; fluid = inst)[:Reservoir]
                @test collect(keys(r.primary_variables)) == [:Pressure, :CO2X, :sGmax]
                @test !haskey(r.equations, :dissolved_co2)
                @test r.equations[:hysteresis].dissolution isa
                      JutulCO2Lab.InstantDissolution
            end
            @test_throws ErrorException co2lab_ve_fluid(default_pvt = nothing,
                                                        co2_rho_pvt = (1, 2, 3, 4))

            # the VE fluid needs column heights
            plain = reservoir_domain(mesh, permeability = 0.5Darcy, porosity = 0.25)
            @test_throws ErrorException setup_co2lab_model(plain)
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

        @testset "hysteresis: sGmax an unknown with co2lab's equation" begin
            hfl = co2lab_ve_fluid(default_pvt = nothing, wat_rho_pvt = (4.3e-10, 100e5),
                                  co2_rho_pvt = (4.3e-9, 100e5), residual = (0.2, 0.2))
            model, states = run(fluid = hfl)
            r = model[:Reservoir]
            @test haskey(r.primary_variables, :sGmax) && !haskey(r.parameters, :sGmax)
            @test r.equations[:hysteresis] isa HysteresisEquation
            @test Jutul.matrix_layout(r.context) isa Jutul.BlockMajorLayout

            # converged sGmax is the running maximum of converged sG
            sG = reduce(hcat, [s[:sG] for s in states])
            sGmax = reduce(hcat, [s[:sGmax] for s in states])
            @test maximum(abs.(sGmax .- accumulate(max, sG, dims = 2))) < 1e-12

            # after shut-in the plume moves on and leaves residual CO2 behind: a
            # residual zone (h < h_max) holding res_gas · (h_max − h) / H of CO2
            hh = states[end][:PlumeHeights]
            zone = findall(i -> hh[2, i] - hh[1, i] > 0.1, axes(hh, 2))
            @test !isempty(zone)
            s_res = [(0.8 * hh[1, i] + 0.2 * (hh[2, i] - hh[1, i])) / 20
                     for i in zone]
            @test sG[zone, end] ≈ s_res
            @test co2_mass(states[end]) ≈ injected rtol = 1e-3
        end

        @testset "vertical-model interface" begin
            # the sharp interface reads no densities; a model through the generic
            # path does, and a stand-in with the same formulas gives the same run
            kw = (wat_rho_pvt = (4.3e-10, 100e5), co2_rho_pvt = (4.3e-9, 100e5),
                  residual = (0.2, 0.2))
            m1, s1 = run(fluid = co2lab_ve_fluid(; default_pvt = nothing, kw...))
            m2, s2 = run(fluid = co2lab_ve_fluid(; default_pvt = nothing,
                                                 vertical = _StandInSharp(), kw...))
            dep(m, k) = Jutul.get_dependencies(m[:Reservoir].secondary_variables[k],
                                               m[:Reservoir])
            for k in (:PlumeHeights, :RelativePermeabilities)
                @test !(:PhaseMassDensities in dep(m1, k))
                @test :PhaseMassDensities in dep(m2, k)
            end
            @test m2[:Reservoir].secondary_variables[:PlumeHeights].vertical isa
                  _StandInSharp
            for k in (:Pressure, :sG, :sGmax, :PlumeHeights)
                @test maximum(abs.(s1[end][k] .- s2[end][k]) ./
                              max.(abs.(s1[end][k]), 1.0)) < 1e-12
            end
        end

        @testset "capillary fringe: runs, conserves CO2" begin
            ff = co2lab_ve_fluid(; default_pvt = nothing,
                                 vertical = CapillaryFringeTable(; invPc3D = (0.4, 0.5),
                                                                   kr3D = 2),
                                 wat_rho_pvt = (4.3e-10, 100e5),
                                 co2_rho_pvt = (4.3e-9, 100e5))
            model, st = run(fluid = ff)
            r = model[:Reservoir]
            @test r.secondary_variables[:PlumeHeights].vertical isa
                  JutulCO2Lab.PScaledTable
            @test abs(co2_mass(st[end]) - injected) / injected < 1e-3
            # the fringe: capillary heights above the sharp interface's for the
            # same saturation, and the capillary pressure Δρ g h
            hh = st[end][:PlumeHeights]
            sharp = JutulCO2Lab._plume_heights.(st[end][:sG], st[end][:sGmax], 20.0,
                                                0.0, 0.0)
            co2 = findall(>(1e-3), st[end][:sG])
            @test !isempty(co2) && all(hh[1, co2] .> first.(sharp[co2]))
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

            # the model's VE variables are the co2lab formulas
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
            dfl = co2lab_ve_fluid(default_pvt = nothing, wat_rho_pvt = (4.3e-10, 100e5),
                                  co2_rho_pvt = (4.3e-9, 100e5), dissolution = true,
                                  dis_rate = 2e-8, dis_max = rs_sat)
            model, st = run(fluid = dfl)
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

    @testset "VE wells" begin
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

        fl = co2lab_ve_fluid(default_pvt = nothing, co2_rho_ref = 700)
        c = co2_injector_control(TotalRateTarget(1.0), fl)
        @test c.injection_mixture == [0.0, 1.0]
        @test c.mixture_density == 700.0
        @test c.phases == ((2, 1.0),)

        f = joinpath(refdir, "sloping_aquifer.mat")
        if !isfile(f)
            @info "Skipping: MRST reference dump $f not found"
        else
            _ve_wells_vs_mrst(matread(f))
        end
    end

    @testset "VE boundary conditions" begin
        f = joinpath(refdir, "ve_bc.mat")
        if !isfile(f)
            @info "Skipping: MRST reference dump $f not found"
        else
            _ve_bc_vs_mrst(matread(f))
        end

        f = joinpath(refdir, "sloping_aquifer.mat")
        if !isfile(f)
            @info "Skipping: MRST reference dump $f not found"
        else
            _ve_bc_runs(matread(f))
        end
    end

    @testset "flow and well model vs MRST runs" begin
        fs = [joinpath(refdir, f) for f in ("sloping_aquifer.mat", "mrst_ref_rate.mat",
                                            "mrst_ref_bhp.mat", "mrst_ref_src.mat")]
        if !all(isfile, fs)
            @info "Skipping: MRST reference runs not found"
        else
            _ve_well_vs_mrst_runs(matread.(fs)...)
        end
    end

    @testset "VE reconstruction" begin
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
            @info "Skipping: MRST reference dump $f not found"
        else
            _ve_reconstruction_vs_mrst(matread(f))
        end
    end

    @testset "capfringeFigures.m, sharp interface, vs MRST" begin
        f = joinpath(refdir, "capfringe_si.mat")
        if !isfile(f)
            @info "Skipping: MRST reference dump $f not found"
        else
            rep = _capfringe_si_vs_mrst(matread(f))
            passed(rep) || display(rep)
            @test passed(rep)
        end
    end

    @testset "capfringeFigures.m, capillary fringe, vs MRST" begin
        f = joinpath(refdir, "capfringe_cp.mat")
        if !isfile(f)
            @info "Skipping: MRST reference dump $f not found"
        else
            rep = _capfringe_cp_vs_mrst(matread(f))
            passed(rep) || display(rep)
            @test passed(rep)
        end
    end

    @testset "capillary fringe ('P-scaled table') vs co2lab" begin
        f = joinpath(refdir, "fringe_fluid.mat")
        if !isfile(f)
            @info "Skipping: MRST reference dump $f not found"
        else
            rep = _fringe_fluid_vs_mrst(matread(f))
            passed(rep) || display(rep)
            @test passed(rep)
        end
    end

    @testset "Leverett-scaled fringe ('P-K-scaled table', corrected)" begin
        J = JutulCO2Lab
        invJ = x -> (0.5 / (x + 0.5))^2
        fl = co2lab_ve_fluid(; default_pvt = nothing,
                             vertical = CapillaryFringeTable(; invPc3D = invJ,
                                                               kr3D = 2,
                                                               leverett = true),
                             co2_rho_ref = 700.0, wat_rho_ref = 1000.0)
        H, drho, g = 20.0, 300.0, Jutul.gravity_constant
        dom = Dict(:column_height => [H, H], :porosity => [0.2, 0.2],
                   :permeability => [1e-13, 4e-13])
        m = J._vertical_model(fl.vertical, dom, fl)
        @test m.kscale ≈ sqrt.(0.2 ./ [1e-13, 4e-13]) .* 30e-3
        # the column's saturation and krG, by direct integration of the fine
        # profile s(z) = 1 − invJ((Pc − Δρ g z)/kscale), thin to truncated plumes
        z = range(0, H; length = 20_001)
        integral(f) = sum(f.(z[1:end-1] .+ step(z) / 2)) * step(z) / H
        for i in 1:2, S in (0.05, 0.3, 0.9)
            h, _ = J._column_heights(m, i, S, S, H, drho, 0.0, 0.0)
            _, krG = J._column_relperm(m, i, S, S, h, h, H, drho, 0.0, 1.0, 1.0, 0.0, 0.0)
            s(zz) = 1 - invJ(max(h - zz, 0.0) * drho * g / m.kscale[i])
            @test integral(s) ≈ S rtol = 1e-4
            @test integral(zz -> s(zz)^2) ≈ krG rtol = 1e-3
        end
    end

    @testset "touched brine of a fringe column (dissolution)" begin
        J = JutulCO2Lab
        H, drho, g, pe = 50.0, 300.0, Jutul.gravity_constant, 5e3
        # with and without entry pressure
        for invpc in (p -> 0.8 * (pe / max(p, pe))^2 + 0.2,
                      p -> 0.8 * (pe / (max(p, 0.0) + pe))^2 + 0.2)
            fl = co2lab_ve_fluid(; default_pvt = nothing,
                                 vertical = CapillaryFringeTable(; invPc3D = invpc),
                                 residual = (0.2, 0.0))
            m = J._vertical_model(fl.vertical, Dict(:column_height => [H]), fl)
            St = m.touched
            hmax_of(s) = J._column_heights(m, 1, s, s, H, drho, 0.2, 0.0)[2]
            touched(s) = 1 - J._untouched_brine(m, 1, s, hmax_of(s), H, drho)
            # above S_t: down to h_max less the entry height
            @test touched(0.1) ≈ (hmax_of(0.1) - m.entry / (drho * g)) / H
            # below S_t: linear in the saturation, continuous at S_t
            @test touched(St / 4) ≈ touched(St) / 4
            @test touched(St * (1 - 1e-9)) ≈ touched(St) rtol = 1e-6
            @test touched(0.0) == 0
        end
        # the sharp interface: the column above h_max, as co2lab
        @test J._untouched_brine(SharpInterface(), 1, 0.1, 5.0, 50.0, 300.0) ≈ 0.9
    end

    @testset "integrated sharp interface" begin
        J = JutulCO2Lab
        # a small sloping aquifer, 4 columns of 10 layers of 2 m
        mesh = _cp_mesh(4, 10, (i, k) -> 1000.0 + 2.0 * i + 2.0 * (k - 1))
        cols, _ = ve_columns(mesh)
        layer = [cell_ijk(mesh, c)[3] for c in 1:number_of_cells(mesh)]
        rw, rg, krw, krg = 0.2, 0.25, 0.75, 0.8

        # uniform rock: the sharp interface's heights and relperm
        u = SharpInterfaceIntegrated(mesh, cols; poro = fill(0.2, 40),
                                     perm_x = fill(1e-13, 40))
        for (s, smax, dh) in ((0.1, 0.1, 0.0), (0.1, 0.4, 0.0), (0.5, 0.6, 1.5))
            hu = J._column_heights(u, 2, s, smax, 20.0, 0.0, rw, rg)
            hs = J._plume_heights(s, smax, 20.0, rw, rg)
            @test all(hu .≈ hs)
            @test all(J._column_relperm(u, 2, s, smax, hu..., 20.0, 0.0, dh, krw, krg,
                                        rw, rg) .≈
                      J._sharp_interface_relperm(hs..., 20.0, krw, krg, dh))
        end

        # layered rock (porosity and permeability by layer)
        poro = [0.1 + 0.02 * l for l in layer]
        perm = [1e-13 * (1 + 0.5 * sin(l)) for l in layer]
        m = SharpInterfaceIntegrated(mesh, cols; poro, perm_x = perm,
                                     perm_y = 2 .* perm)
        @test_logs (:warn, r"integrated relperm") SharpInterfaceIntegrated(mesh, cols;
            poro, perm_x = perm, perm_y = 2 .* perm)
        # the column's fine profile: layer of depth z, and its φ and k
        lay(zz) = clamp(ceil(Int, zz / 2), 1, 10)
        ph(zz) = 0.1 + 0.02 * lay(zz)
        kx(zz) = 1e-13 * (1 + 0.5 * sin(lay(zz)))
        z = range(0, 20.0; length = 200_001)
        mid = z[1:end-1] .+ step(z) / 2
        quad(f, a, b) = sum((f(x) for x in mid if a <= x < b); init = 0.0) * step(z)
        # states with sG at least the residual left by sGmax, sGmax·rg/(1 − rw)
        for (s, smax, dh) in ((0.05, 0.05, 0.0), (0.2, 0.4, 0.0), (0.5, 0.6, 1.5))
            h, hmax = J._column_heights(m, 2, s, smax, 20.0, 0.0, rw, rg)
            # pore volumes: CO2 at 1 − rw in the plume, rg in the residual zone
            Vt = quad(ph, 0, 20)
            @test (1 - rw) * quad(ph, 0, hmax) ≈ smax * Vt rtol = 1e-4
            @test (1 - rw - rg) * quad(ph, 0, h) + rg * quad(ph, 0, hmax) ≈
                  s * Vt rtol = 1e-4
            # relperm by direct integration of the zones
            krW, krG = J._column_relperm(m, 2, s, smax, h, hmax, 20.0, 0.0, dh, krw,
                                         krg, rw, rg)
            KH = quad(kx, 0, 20)
            @test krG ≈ krg * quad(kx, min(dh, h), h) / KH rtol = 1e-4
            @test krW ≈ (quad(kx, max(hmax, dh), 20) +
                         krw * quad(kx, max(h, dh), max(hmax, dh))) / KH rtol = 1e-4
            # the touched brine: pore volume above h_max
            @test 1 - J._untouched_brine(m, 2, smax, hmax, 20.0, 0.0) ≈
                  quad(ph, 0, hmax) / Vt rtol = 1e-4
        end

        # the inventory split: co2lab's by height; ours by pore volume, whose
        # parts add up as the height split's do
        S, Smax = fill(0.3, 4), fill(0.45, 4)
        hh = [J._column_heights(m, i, 0.3, 0.45, 20.0, 0.0, rw, rg) for i in 1:4]
        col = (H = fill(20.0, 4), zt = fill(5.0, 4), dh = fill(1.0, 4))
        split(c) = J._category_saturations(m, S, Smax, first.(hh), last.(hh), col,
                                           zeros(4), rw, rg; co2lab_compatible = c)
        p, c = split(false), split(true)
        @test c == J._height_split(S, Smax, first.(hh), last.(hh), col, rw, rg)
        S_immob = Smax .* rg ./ (1 - rw)
        @test p.immob_sub .+ p.immob_struct .+ p.immob_free ≈ S_immob
        @test p.mobile_sub .+ p.mobile_struct .+ p.mobile_free ≈ S .- S_immob
        # porosity grows with depth, so the top 5 m hold less than their height share
        @test all(p.mobile_struct .< c.mobile_struct)

        # a run: CO2 conserved, and no dependency on the densities
        ve = ve_domain(mesh, cols; perm_x = perm, poro)
        fl = co2lab_ve_fluid(; default_pvt = nothing, vertical = m, residual = (rw, rg),
                             wat_rho_pvt = (4.3e-10, 100e5),
                             co2_rho_pvt = (4.3e-9, 100e5))
        well = setup_ve_well(ve, 4; name = :inj)
        model = setup_co2lab_model(ve.domain; wells = [well], fluid = fl)
        q = 2e-3
        on = setup_reservoir_forces(model;
                                    control = Dict(:inj => co2_injector_control(
                                        TotalRateTarget(q), fl)))
        off = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()))
        p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = 1100.0)
        state0 = setup_co2lab_state(model; pressure = p0, forces = on)
        dt = fill(30 * 86400.0, 10)
        case = JutulCase(model, dt, [k <= 5 ? on : off for k in 1:10]; state0,
                         parameters = setup_parameters(model))
        sim, cfg = setup_co2lab_simulator(case; info_level = -1)
        st = simulate_reservoir(case; simulator = sim, config = cfg).states
        injected = 5 * 30 * 86400.0 * q * fl.rho_ref[2]
        # (the sharp interface misses by the same 1.8e-5 on this case)
        @test abs(sum(st[end][:TotalMasses][2, :]) - injected) / injected < 1e-4
        r = model[:Reservoir]
        @test !(:PhaseMassDensities in
                Jutul.get_dependencies(r.secondary_variables[:PlumeHeights], r))
        # a model built for another grid is refused
        small = _cp_mesh(2, 10, (i, k) -> 1000.0 + 2.0 * (k - 1))
        other = SharpInterfaceIntegrated(small, ve_columns(small)[1];
                                         poro = fill(0.2, 20), perm_x = fill(1e-13, 20))
        @test_throws ErrorException setup_co2lab_model(ve.domain;
            fluid = co2lab_ve_fluid(; default_pvt = nothing, vertical = other))
    end

    @testset "layered fringe" begin
        J = JutulCO2Lab
        # 4 columns of 10 layers of 2 m; Δρ, residuals; fringes with and
        # without entry pressure
        mesh = _cp_mesh(4, 10, (i, k) -> 1000.0 + 2.0 * i + 2.0 * (k - 1))
        cols, _ = ve_columns(mesh)
        layer = [cell_ijk(mesh, c)[3] for c in 1:number_of_cells(mesh)]
        H, drho, g, rw, rg, pe = 20.0, 300.0, Jutul.gravity_constant, 0.2, 0.25, 3e3
        curves = (p -> 0.8 * (pe / (max(p, 0.0) + pe))^2 + 0.2,
                  p -> 0.8 * (pe / max(p, pe))^2 + 0.2)
        dom = Dict(:column_height => fill(H, 4))
        fl(v) = co2lab_ve_fluid(; default_pvt = nothing, vertical = v, residual = (rw, rg))
        build(v) = J._vertical_model(v, dom, fl(v))
        uniform = (poro = fill(0.2, 40), perm_x = fill(1e-13, 40))

        # uniform rock: CapillaryFringeTable's heights and krG above the tip;
        # krW co2lab's krW_simple less r_w/(1 − r_w)·(S − S_free) (the
        # integrated brine)
        for (invpc, entry) in zip(curves, (0.0, pe))
            mu = build(LayeredFringe(mesh, cols; uniform..., invPc3D = invpc, kr3D = 2))
            mc = build(CapillaryFringeTable(; invPc3D = invpc, kr3D = 2))

            # the tip: tip_height·H below the entry height; below it the height
            # linear in the saturation, meeting the exact one at S_t
            ht, St = J._tip(mu, 2, drho)
            @test ht ≈ entry / (drho * g) + 1e-6 * H rtol = 1e-12
            @test 0 < St < 1e-5
            @test J._layered_height(mu, 2, St / 3, drho) ≈ ht / 3 rtol = 1e-12
            @test J._layered_newton(mu, 2, St, drho) ≈ ht rtol = 1e-8

            for (s, smax) in ((0.05, 0.05), (0.2, 0.4), (0.79, 0.79))
                hu = J._column_heights(mu, 2, s, smax, H, drho, rw, rg)
                hc = J._column_heights(mc, 2, s, smax, H, drho, rw, rg)
                @test all(isapprox.(hu, hc; rtol = 1e-12))
                ku = J._column_relperm(mu, 2, s, smax, hu..., H, drho, 0.0, 1.0, 1.0,
                                       rw, rg)
                kc = J._column_relperm(mc, 2, s, smax, hc..., H, drho, 0.0, 1.0, 1.0,
                                       rw, rg)
                sf = J._free_saturation(s, smax, rw, rg)
                @test ku[2] ≈ kc[2] rtol = 1e-12
                @test ku[1] ≈ kc[1] - rw / (1 - rw) * (s - sf) atol = 1e-12
            end
        end

        # layered rock: saturations by pore volume, relperms by permeability,
        # against quadrature of the fine profile; an interval inside the column
        poro = [0.1 + 0.02 * l for l in layer]
        perm = [1e-13 * (1 + 0.5 * sin(l)) for l in layer]
        lay(zz) = clamp(ceil(Int, zz / 2), 1, 10)
        ph(zz) = 0.1 + 0.02 * lay(zz)
        kx(zz) = 1e-13 * (1 + 0.5 * sin(lay(zz)))
        z = range(0, H; length = 200_001)
        mid = z[1:end-1] .+ step(z) / 2
        quad(f, a, b) = sum((f(x) for x in mid if a <= x < b); init = 0.0) * step(z)
        C = rg / (1 - rw)
        quadratic = sw -> ((sw - rw) / (1 - rw))^2
        for invpc in curves, krw3D in (:linear, quadratic)
            m = build(LayeredFringe(mesh, cols; poro, perm_x = perm, invPc3D = invpc,
                                    kr3D = 2, krw3D))
            kw = krw3D isa Symbol ? (sw -> (sw - rw) / (1 - rw)) : krw3D
            for (s, smax) in ((0.05, 0.05), (0.2, 0.4), (0.5, 0.6))
                h, hmax = J._column_heights(m, 2, s, smax, H, drho, rw, rg)
                krW, krG = J._column_relperm(m, 2, s, smax, h, hmax, H, drho, 0.0, 1.0,
                                             1.0, rw, rg)
                sg(p) = 1 - invpc(max(p, 0.0))
                sf(zz) = sg(drho * g * (h - zz))
                sa(zz) = (1 - C) * sf(zz) + C * sg(drho * g * (hmax - zz))
                Vt, KH = quad(ph, 0, H), quad(kx, 0, H)
                @test quad(zz -> ph(zz) * sg(drho * g * (hmax - zz)), 0, H) / Vt ≈
                      smax rtol = 1e-6
                @test quad(zz -> ph(zz) * sa(zz), 0, H) / Vt ≈ s rtol = 1e-6
                @test quad(zz -> kx(zz) * sf(zz)^2, 0, H) / KH ≈ krG rtol = 1e-5
                @test quad(zz -> kx(zz) * kw(1 - sa(zz)), 0, H) / KH ≈ krW rtol = 1e-6
                # the cells are 100 m × 100 m in plan (relative to the column's)
                I = J._interval_integrals(m, 2, 3.3, 11.7, h, hmax, drho)
                @test I.co2 ≈ 1e4 * quad(zz -> ph(zz) * sa(zz), 3.3, 11.7) atol =
                      1e-6 * s * 1e4 * Vt
            end
        end

        # zero capillarity: the integrated sharp interface (its brine endpoint
        # in the residual zone that of the linear krw3D, 1 − r_g/(1 − r_w))
        sharp = build(LayeredFringe(mesh, cols; poro, perm_x = perm,
                                    invPc3D = p -> p > 0 ? rw : 1.0, kr3D = 2))
        si = SharpInterfaceIntegrated(mesh, cols; poro, perm_x = perm)
        for (s, smax) in ((0.05, 0.05), (0.2, 0.4))
            hl = J._column_heights(sharp, 2, s, smax, H, drho, rw, rg)
            hs = J._column_heights(si, 2, s, smax, H, drho, rw, rg)
            @test all(isapprox.(hl, hs; rtol = 1e-9))
            kl = J._column_relperm(sharp, 2, s, smax, hl..., H, drho, 0.0, 1.0, 1.0,
                                   rw, rg)
            ks = J._column_relperm(si, 2, s, smax, hs..., H, drho, 0.0,
                                   1 - rg / (1 - rw), (1 - rw)^2, rw, rg)
            @test all(isapprox.(kl, ks; rtol = 1e-9))
        end
        # ... with its derivatives, also at an empty column's S = 0 (the
        # profile steps up at p = 0; a quadratic brine curve integrated by
        # quadrature across the moving step)
        sq = build(LayeredFringe(mesh, cols; poro, perm_x = perm,
                                 invPc3D = p -> p > 0 ? rw : 1.0, kr3D = 2,
                                 krw3D = sw -> ((sw - rw) / (1 - rw))^2))
        # no residual CO2 (rg = 0): the brine endpoints of both are 1
        kr(m, x, kw) = (sm = J._effective_max(x[2], x[1], rw, 0.0);
                        hh = J._column_heights(m, 2, x[1], sm, H, drho, rw, 0.0);
                        collect(J._column_relperm(m, 2, x[1], sm, hh..., H, drho, 0.0,
                                                  kw, (1 - rw)^2, rw, 0.0)))
        for x0 in ([0.0, 0.0], [0.1, 0.0])
            jl = Jutul.ForwardDiff.jacobian(x -> kr(sq, x, 1.0), x0)
            js = Jutul.ForwardDiff.jacobian(x -> kr(si, x, 1.0), x0)
            @test jl ≈ js rtol = 1e-8
        end

        # derivatives (one Newton step with AD) against central differences,
        # away from the drainage–imbibition switch
        m = build(LayeredFringe(mesh, cols; poro, perm_x = perm, invPc3D = curves[2],
                                kr3D = 2))
        f(x) = (hh = J._column_heights(m, 3, x[1], x[2], H, x[3], rw, rg);
                [hh..., J._column_relperm(m, 3, x[1], x[2], hh..., H, x[3], 0.0, 1.0,
                                          1.0, rw, rg)...])
        for x0 in ([0.2, 0.4, 300.0], [0.04, 0.05, 280.0])
            jac = Jutul.ForwardDiff.jacobian(f, x0)
            e = [1e-7, 1e-7, 1e-4] .* x0
            fd = reduce(hcat, [(f(x0 .+ e[j] .* (1:3 .== j)) .-
                                f(x0 .- e[j] .* (1:3 .== j))) ./ (2e[j]) for j in 1:3])
            @test jac ≈ fd rtol = 1e-5
        end
        # an empty column of uniform rock keeps co2lab's brine slope,
        # −1/(1 − r_w), at S = 0 (layers weight it by the top's k/φ)
        mu = build(LayeredFringe(mesh, cols; uniform..., invPc3D = curves[2], kr3D = 2))
        @test Jutul.ForwardDiff.derivative(
            s -> J._column_relperm(mu, 3, s, s, J._column_heights(mu, 3, s, s, H, drho,
                                   rw, rg)..., H, drho, 0.0, 1.0, 1.0, rw, rg)[1],
            0.0) ≈ -1 / (1 - rw) rtol = 1e-6

        # co2lab's dimensionless curve, physical: CapillaryFringeTable's scaling
        cc = co2lab_capillary_curve(0.4, 0.5; height = H, drho = 300.0, res_water = rw)
        scaled = J._fine_capillary_curve((0.4, 0.5), H * 300.0 * g, rw)
        @test all(cc(p) ≈ scaled(p) for p in (0.0, 1e3, 1e5))
        @test_throws ArgumentError LayeredFringe(mesh, cols; poro, perm_x = perm,
                                                 invPc3D = (0.4, 0.5))
        # co2lab's fringe refuses the caprock's rugosity; LayeredFringe takes it
        rough = Dict(:column_height => fill(H, 4), :rugosity => fill(0.5, 4))
        v = CapillaryFringeTable(; invPc3D = curves[1])
        @test_throws ArgumentError J._vertical_model(v, rough, fl(v))
        v = LayeredFringe(mesh, cols; poro, perm_x = perm, invPc3D = curves[1])
        @test J._vertical_model(v, rough, fl(v)) isa J.LayeredFringeModel

        # rugosity: the relief's CO2 immobile, both relperms integrated from dh
        # (quadrature), at zero capillarity the integrated sharp interface's
        # accretion layer, with its derivatives
        for invpc in curves, (s, smax, dh) in ((0.2, 0.4, 1.5), (0.5, 0.6, 3.0))
            m = build(LayeredFringe(mesh, cols; poro, perm_x = perm, invPc3D = invpc,
                                    kr3D = 2))
            h, hmax = J._column_heights(m, 2, s, smax, H, drho, rw, rg)
            krW, krG = J._column_relperm(m, 2, s, smax, h, hmax, H, drho, dh, 1.0, 1.0,
                                         rw, rg)
            sg(p) = 1 - invpc(max(p, 0.0))
            sf(zz) = sg(drho * g * (h - zz))
            sa(zz) = (1 - C) * sf(zz) + C * sg(drho * g * (hmax - zz))
            KH = quad(kx, 0, H)
            @test quad(zz -> kx(zz) * sf(zz)^2, dh, H) / KH ≈ krG rtol = 1e-5
            @test quad(zz -> kx(zz) * (1 - sa(zz) / (1 - rw)), dh, H) / KH ≈ krW rtol = 1e-6
        end
        rel(m, kw, d) = x -> (sm = J._effective_max(x[2], x[1], rw, rg);
                              hh = J._column_heights(m, 2, x[1], sm, H, drho, rw, rg);
                              collect(J._column_relperm(m, 2, x[1], sm, hh..., H, drho,
                                                        d, kw, (1 - rw)^2, rw, rg)))
        # (off the drainage–imbibition switch, where the two models take
        # different one-sided derivatives)
        for (s, smax, d) in ((0.04, 0.05, 1.5), (0.2, 0.4, 3.0), (0.5, 0.6, 0.5))
            x0, kwi = [s, smax], 1 - rg / (1 - rw)
            @test rel(sharp, 1.0, d)(x0) ≈ rel(si, kwi, d)(x0) rtol = 1e-9
            @test Jutul.ForwardDiff.jacobian(rel(sharp, 1.0, d), x0) ≈
                  Jutul.ForwardDiff.jacobian(rel(si, kwi, d), x0) rtol = 1e-8
        end

        # a run, with rugosity: CO2 conserved
        ve = ve_domain(mesh, cols; perm_x = perm, poro, rugosity = 0.5)
        v = LayeredFringe(mesh, cols; poro, perm_x = perm, invPc3D = curves[1], kr3D = 2,
                          krw3D = quadratic)
        flr = co2lab_ve_fluid(; default_pvt = nothing, vertical = v, residual = (rw, rg),
                              wat_rho_pvt = (4.3e-10, 100e5),
                              co2_rho_pvt = (4.3e-9, 100e5))
        well = setup_ve_well(ve, 4; name = :inj)
        model = setup_co2lab_model(ve.domain; wells = [well], fluid = flr)
        q = 2e-4
        on = setup_reservoir_forces(model;
                                    control = Dict(:inj => co2_injector_control(
                                        TotalRateTarget(q), flr)))
        off = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()))
        p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = 1100.0)
        state0 = setup_co2lab_state(model; pressure = p0, forces = on)
        dt = fill(30 * 86400.0, 10)
        case = JutulCase(model, dt, [k <= 5 ? on : off for k in 1:10]; state0,
                         parameters = setup_parameters(model))
        sim, cfg = setup_co2lab_simulator(case; info_level = -1)
        st = simulate_reservoir(case; simulator = sim, config = cfg).states
        injected = 5 * 30 * 86400.0 * q * flr.rho_ref[2]
        @test abs(sum(st[end][:TotalMasses][2, :]) - injected) / injected < 1e-4
        @test model[:Reservoir].secondary_variables[:PlumeHeights].vertical isa
              J.LayeredFringeModel
    end

    @testset "fine-scale fluid and fine cells" begin
        J = JutulCO2Lab
        mesh = _cp_mesh_pillars(4, 6, (p, k) -> 1000.0 - 2.0 * p + 1.0 * (k - 1))
        nc = number_of_cells(mesh)
        layer = [cell_ijk(mesh, c)[3] for c in 1:nc]
        perm = [1e-13 * (1 + 0.5 * sin(l)) for l in layer]
        poro = [0.15 + 0.01 * l for l in layer]
        rw, rg, pe, g = 0.2, 0.2, 2e3, Jutul.gravity_constant
        sand = p -> 0.8 * (pe / (max(p, 0.0) + pe))^2 + 0.2
        # a tighter rock with an entry pressure, in layers 3 and 4
        shale = p -> 0.8 * (4pe / max(p, 4pe))^2 + 0.2
        krg(s) = min(s, 1 - rw)^2
        krw(sw) = (max(sw - rw, 0.0) / (1 - rw))^2
        rt = [l in (3, 4) ? 2 : 1 for l in layer]
        ff = FineScaleFluid(; rock_types = [RockCurves(; invPc3D = sand, kr3D = krg,
                                                         krw3D = krw),
                                            RockCurves(; invPc3D = shale, kr3D = 3)],
                            residual = (rw, rg), default_pvt = nothing,
                            co2_rho_pvt = nothing, wat_rho_pvt = (4.3e-10, 1e7))
        cols, _ = ve_columns(mesh)
        fl = co2lab_ve_fluid(ff, mesh, cols; poro, perm_x = perm, rock_type = rt)
        @test fl.vertical isa LayeredFringe && fl.fine.rock_type == rt

        # the VE column of two rock types against quadrature (1 m layers)
        dom = Dict(:column_height => fill(6.0, 4))
        m = J._vertical_model(fl.vertical, dom, fl)
        drho = 300.0
        h, hmax = J._column_heights(m, 2, 0.2, 0.3, 6.0, drho, rw, rg)
        typ(z) = clamp(ceil(Int, z), 1, 6) in (3, 4) ? 2 : 1
        ph(z) = 0.15 + 0.01 * clamp(ceil(Int, z), 1, 6)
        sg(z, hh) = 1 - (typ(z) == 1 ? sand : shale)(max(drho * g * (hh - z), 0.0))
        C = rg / (1 - rw)
        zs = range(0, 6; length = 600_001)
        mid = zs[1:end-1] .+ step(zs) / 2
        quad(f) = sum(f(z) for z in mid) * step(zs)
        @test quad(z -> ph(z) * ((1 - C) * sg(z, h) + C * sg(z, hmax))) / quad(ph) ≈ 0.2 rtol = 1e-5
        @test quad(z -> ph(z) * sg(z, hmax)) / quad(ph) ≈ 0.3 rtol = 1e-5

        # without capillarity: the integrated sharp interface, the curves'
        # endpoints; rock types with different endpoints refused
        flat = FineScaleFluid(; rock_types = [RockCurves(; kr3D = krg, krw3D = krw)],
                              residual = (rw, rg), default_pvt = nothing)
        fs = co2lab_ve_fluid(flat, mesh, cols; poro, perm_x = perm)
        @test fs.vertical isa SharpInterfaceIntegrated
        @test fs.krg ≈ krg(1 - rw) && fs.krw ≈ krw(1 - rg)
        two = FineScaleFluid(; rock_types = [RockCurves(kr3D = 2), RockCurves(kr3D = 3)],
                             residual = (rw, rg), default_pvt = nothing)
        @test_throws ArgumentError co2lab_ve_fluid(two, mesh, cols; poro, perm_x = perm,
                                                   rock_type = rt)

        # fine cells: relperm at the free saturation, the brine at 1 − sG; the
        # capillary pressure at the free saturation, the entry pressure at none
        hcols, _ = ve_columns(mesh; is_ve = fill(false, nc), below_holes = :stack)
        hve = ve_domain(mesh, hcols; perm_x = perm, poro, fault_flux = :virtual)
        hfl = co2lab_ve_fluid(ff, mesh, hcols; poro, perm_x = perm, rock_type = rt)
        f = J._fine_cell_fluid(hfl.fine, hve.domain)
        b2 = findfirst(==(2), f.rock_type)
        b1 = findfirst(==(1), f.rock_type)
        sf = J._free_saturation(0.3, 0.5, rw, rg)
        @test all(isapprox.(J._fine_relperm(f, b1, 0.3, 0.5), (krw(0.7), krg(sf));
                            atol = 1e-6))
        @test J._fine_capillary_pressure(f, b2, 0.0, 0.0) ≈ 4pe rtol = 1e-3
        p1 = J._fine_capillary_pressure(f, b1, 0.3, 0.3)
        @test 1 - sand(p1) ≈ 0.3 rtol = 1e-3

        # an all-fine grid runs, and reproduces JutulDarcy's 3D model with the
        # same curves (to the formulations' difference: linear against
        # exponential brine compressibility, black-oil against immiscible)
        RHOW, RHOC = 1050.0, 700.0
        ff1 = FineScaleFluid(; rock_types = [RockCurves(; invPc3D = sand, kr3D = krg,
                                                          krw3D = krw)],
                             residual = (rw, 0.0), default_pvt = nothing,
                             co2_rho_pvt = nothing, co2_rho_ref = RHOC,
                             wat_rho_pvt = (4.3e-10, 1e7), wat_rho_ref = RHOW,
                             co2_mu_pvt = nothing, co2_mu_ref = 6e-5,
                             wat_mu_pvt = nothing, wat_mu_ref = 8e-4)
        hf1 = co2lab_ve_fluid(ff1, mesh, hcols; poro, perm_x = perm)
        cell(i, k) = findfirst(c -> cell_ijk(mesh, c) == (i, 1, k), 1:nc)
        ic, pc_ = cell(1, 6), cell(4, 1)
        blk(c) = hcols.partition[c]
        z = hve.domain[:cell_centroids][3, :]
        p0 = RHOW * g .* z
        bhp, q = p0[blk(pc_)], 5e-5
        dt = fill(2 * 86400.0, 8)
        kw = (info_level = -1, tol_cnv = 1e-10, tol_mb = 1e-12, timesteps = :none,
              max_timestep = Inf)
        wells = [setup_well(hve.domain, [blk(ic)]; name = :inj, WI = 1e-12),
                 setup_well(hve.domain, [blk(pc_)]; name = :prod, WI = 1e-12)]
        model = setup_co2lab_model(hve.domain; wells, fluid = hf1,
                                   rock_compressibility = nothing)
        ctrl = Dict(:inj => co2_injector_control(TotalRateTarget(q), hf1),
                    :prod => ProducerControl(BottomHolePressureTarget(bhp)))
        forces = setup_reservoir_forces(model; control = ctrl)
        state0 = setup_co2lab_state(model; pressure = p0, forces)
        case = JutulCase(model, dt, forces; state0, parameters = setup_parameters(model))
        sim, cfg = setup_co2lab_simulator(case; kw...)
        sh = simulate_reservoir(case; simulator = sim, config = cfg).states

        fine = reservoir_domain(mesh; permeability = perm, porosity = poro)
        sys = ImmiscibleSystem((AqueousPhase(), VaporPhase());
                               reference_densities = (RHOW, RHOC))
        m3 = setup_reservoir_model(fine, sys;
                                   wells = [setup_well(fine, [ic]; name = :inj, WI = 1e-12),
                                            setup_well(fine, [pc_]; name = :prod, WI = 1e-12)])
        s = collect(range(0.0, 1.0; length = 2001))
        replace_variables!(m3; RelativePermeabilities = JutulDarcy.ReservoirRelativePermeabilities(
                               w = JutulDarcy.PhaseRelativePermeability(s, krw.(s); label = :w),
                               g = JutulDarcy.PhaseRelativePermeability(s, krg.(s); label = :g)),
                           PhaseMassDensities = ConstantCompressibilityDensities(
                               p_ref = 1e7, density_ref = [RHOW, RHOC],
                               compressibility = [4.3e-10, 0.0]))
        pcf = J._capillary_of_co2(sand, rw)
        set_secondary_variables!(m3[:Reservoir]; CapillaryPressure =
            JutulDarcy.SimpleCapillaryPressure(Jutul.get_1d_interpolator(s, pcf.(s))))
        f3 = setup_reservoir_forces(m3; control = Dict(
            :inj => InjectorControl(TotalRateTarget(q), [0.0, 1.0], density = RHOC),
            :prod => ProducerControl(BottomHolePressureTarget(bhp))))
        prm = setup_parameters(m3)
        prm[:Reservoir][:PhaseViscosities][1, :] .= 8e-4
        prm[:Reservoir][:PhaseViscosities][2, :] .= 6e-5
        s3 = simulate_reservoir(setup_reservoir_state(m3; Pressure = RHOW * g .*
                                                      fine[:cell_centroids][3, :],
                                                      Saturations = [1.0, 0.0]),
                                m3, dt; parameters = prm, forces = f3, kw...).states
        cellof = [hcols.cells[hcols.pos[b]] for b in 1:number_of_columns(hcols)]
        @test maximum(sh[end][:sG]) > 0.01
        @test maximum(abs.(sh[end][:sG] .- s3[end][:Saturations][2, cellof])) < 1e-4
        @test maximum(abs.(sh[end][:Pressure] .- s3[end][:Pressure][cellof])) <
              1e-5 * maximum(s3[end][:Pressure])

        # instantaneous dissolution in fine cells: the brine saturated wherever
        # CO2 is present, under co2lab's rate law for VE blocks
        ffd = FineScaleFluid(; rock_types = [RockCurves(; invPc3D = sand, kr3D = krg)],
                             residual = (rw, 0.0), default_pvt = nothing,
                             co2_rho_pvt = nothing, wat_rho_pvt = (4.3e-10, 1e7),
                             dissolution = true, dis_max = 0.03)
        hfd = co2lab_ve_fluid(ffd, mesh, hcols; poro, perm_x = perm, dis_rate = 1e-12)
        md = setup_co2lab_model(hve.domain; wells, fluid = hfd)
        ctrl = Dict(:inj => co2_injector_control(TotalRateTarget(q), hfd),
                    :prod => ProducerControl(BottomHolePressureTarget(bhp)))
        fd = setup_reservoir_forces(md; control = ctrl)
        cd = JutulCase(md, dt, fd; state0 = setup_co2lab_state(md; pressure = p0,
                                                                forces = fd),
                       parameters = setup_parameters(md))
        sim, cfg = setup_co2lab_simulator(cd; info_level = -1)
        sd = simulate_reservoir(cd; simulator = sim, config = cfg).states[end]
        co2 = findall(>(1e-3), sd[:sG])
        @test !isempty(co2) && all(abs.(sd[:Rs][co2] .- 0.03) .< 1e-6)
    end

    @testset "transition flux" begin
        J = JutulCO2Lab
        g = Jutul.gravity_constant
        rw, RHOW, RHOC = 0.2, 1050.0, 700.0
        drho = RHOW - RHOC
        # nearly incompressible brine by default (its hydrostatic pressure
        # linear, for the equilibrium checks); closed domains with injection
        # need a compressible one
        fine_fluid(curve; cw = 4.3e-14) =
            FineScaleFluid(; rock_types = [RockCurves(; invPc3D = curve, kr3D = 2)],
                           residual = (rw, 0.0), default_pvt = nothing,
                           co2_rho_pvt = nothing, co2_rho_ref = RHOC,
                           wat_rho_pvt = (cw, 1e7), wat_rho_ref = RHOW)
        # the CO2 saturation of fine cell b at capillary pressure pc (bisection
        # on the cell's own curve; zero below the entry pressure)
        function fine_sg(f, b, pc)
            pcf(s) = J._fine_capillary_pressure(f, b, s, s)
            pcf(0.0) >= pc && return 0.0
            lo, hi = 0.0, 1 - rw
            for _ in 1:100
                mid = (lo + hi) / 2
                pcf(mid) < pc ? (lo = mid) : (hi = mid)
            end
            return (lo + hi) / 2
        end

        # P1: a VE column beside a fine column, both at capillary–gravity
        # equilibrium with a 3 m plume, stays at rest (zero and nonzero entry
        # pressure)
        mesh = _cp_mesh_pillars(2, 12, (p, k) -> 1000.0 + 0.5 * (k - 1))
        nc = number_of_cells(mesh)
        poro, perm = fill(0.2, nc), fill(1e-13, nc)
        cols, _ = ve_columns(mesh; is_ve = [cell_ijk(mesh, c)[1] == 2 for c in 1:nc])
        ve = ve_domain(mesh, cols; perm_x = perm, poro, fault_flux = :virtual)
        @test count(==(J._FACE_VE_FINE), ve.domain[:face_category]) == 12
        P, h = 1e7, 3.0
        for curve in (p -> 0.8 * (2e3 / (max(p, 0.0) + 2e3))^2 + 0.2,
                      p -> 0.8 * (1e3 / max(p, 1e3))^2 + 0.2)
            fl = co2lab_ve_fluid(fine_fluid(curve), mesh, cols; poro, perm_x = perm)
            model = setup_co2lab_model(ve.domain; fluid = fl,
                                       rock_compressibility = nothing)
            f = J._fine_cell_fluid(fl.fine, ve.domain)
            vb = findfirst(cols.is_ve)
            m = model[:Reservoir].secondary_variables[:PlumeHeights].vertical
            H = J._column_extent(m.cells, vb)[2]
            x = ve.domain[:cell_centroids][3, :] .- ve.domain[:reference_depth][vb]
            sg = [b == vb ? J._interval_integrals(m, vb, 0.0, H, h, h, drho).free / H :
                  fine_sg(f, b, drho * g * (h - x[b])) for b in eachindex(x)]
            p = [b == vb ? P : P + RHOW * g * x[b] for b in eachindex(x)]
            state0 = setup_co2lab_state(model; pressure = p, sG = sg)
            case = JutulCase(model, fill(30 * 86400.0, 3), setup_reservoir_forces(model);
                             state0, parameters = setup_parameters(model))
            sim, cfg = setup_co2lab_simulator(case; info_level = -1, tol_cnv = 1e-10,
                                              tol_mb = 1e-10)
            st = simulate_reservoir(case; simulator = sim, config = cfg).states
            @test maximum(sg) > 0.3
            @test maximum(abs.(st[end][:sG] .- sg)) < 1e-8
            # reconstruction: the fine cells keep their own state
            rec = ve_reconstruct(st[end], ve, cols, fl)
            fb = findall(.!cols.is_ve)
            fc = cols.cells[cols.pos[fb]]
            @test rec.sG[fc] == st[end][:sG][fb]
            @test rec.pressure[fc] == st[end][:Pressure][fb]
            @test maximum(abs.(st[end][:Pressure] .- p)) < 1e-3
        end

        # open hydrostatic boundaries on the fine column's side (from each
        # cell's centroid) and the VE column's: brine at rest stays at rest
        sides, _ = ve_boundary_sides(mesh, cols, ve.domain)
        ix = ve_select_sides(sides, cols; directions = (:left, :right))
        @test length(ix) == 13
        fl = co2lab_ve_fluid(fine_fluid(p -> 0.8 * (2e3 / (max(p, 0.0) + 2e3))^2 + 0.2),
                             mesh, cols; poro, perm_x = perm)
        model = setup_co2lab_model(ve.domain; fluid = fl, rock_compressibility = nothing)
        hyd(z) = P + RHOW * g * (z - 1000.0)
        p = hyd.(ve.domain[:reference_depth])
        forces = setup_reservoir_forces(model;
                                        bc = ve_boundary_conditions(ve, sides, ix;
                                                                    pressure = hyd))
        case = JutulCase(model, fill(30 * 86400.0, 3), forces;
                         state0 = setup_co2lab_state(model; pressure = p),
                         parameters = setup_parameters(model))
        sim, cfg = setup_co2lab_simulator(case; info_level = -1)
        st = simulate_reservoir(case; simulator = sim, config = cfg).states
        @test maximum(abs.(st[end][:Pressure] .- p)) < 1e-3

        # a thin plume under a fault offset larger than itself: co2lab's smooth
        # flux spreads it across, the virtual cell keeps it (the caprock up-dip
        # is below the plume); CO2 conserved either way
        mesh = _cp_mesh(8, 6, (i, k) -> 1000.0 - 1.0 * i + (i > 4 ? 2.5 : 0.0) +
                                         1.0 * (k - 1))
        nc = number_of_cells(mesh)
        cols, _ = ve_columns(mesh)
        fl = co2lab_ve_fluid(; default_pvt = nothing, residual = (0.2, 0.2),
                             wat_rho_pvt = (4.3e-10, 1e7), co2_rho_pvt = (4.3e-9, 1e7))
        q, nt = 2e-5, 20
        sG = map((:co2lab, :virtual)) do ff
            ve = ve_domain(mesh, cols; perm_x = fill(1e-13, nc), poro = fill(0.2, nc),
                           fault_flux = ff)
            model = setup_co2lab_model(ve.domain; fluid = fl,
                                       wells = [setup_ve_well(ve, 8; name = :inj)])
            on = setup_reservoir_forces(model; control = Dict(
                :inj => co2_injector_control(TotalRateTarget(q), fl)))
            off = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()))
            p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = 1100.0)
            case = JutulCase(model, fill(30 * 86400.0, nt),
                             [k <= nt ÷ 2 ? on : off for k in 1:nt];
                             state0 = setup_co2lab_state(model; pressure = p0, forces = on),
                             parameters = setup_parameters(model))
            sim, cfg = setup_co2lab_simulator(case; info_level = -1)
            st = simulate_reservoir(case; simulator = sim, config = cfg).states
            injected = nt ÷ 2 * 30 * 86400.0 * q * fl.rho_ref[2]
            @test abs(sum(st[end][:TotalMasses][2, :]) - injected) / injected < 1e-4
            st[end][:sG]
        end
        @test maximum(sG[1][1:6]) > 1e-4
        @test maximum(sG[2][1:6]) < 1e-6 && sG[2][8] > maximum(sG[1][7:8])

        # a fine column inside a sloping VE aquifer (20 layers of 0.3 m): CO2
        # injected down-dip migrates through it into the VE columns up-dip.
        # Under VE conditions fine cells approach VE only as they resolve the
        # plume (0.6 m fringe): the hybrid lies between all-fine and all-VE and
        # converges to VE with refinement (CO2 up-dip after 40 months, 6 / 24 /
        # 60 / 120 layers: 0.17 / 0.76 / 1.8 / 2.25 t, VE 2.39 t; all-fine
        # 2e-5 / 0.23 / 1.5 / 2.1 t). CO2 conserved.
        mesh = _cp_mesh_pillars(6, 20, (p, k) -> 1000.0 + 2.0 * p + 0.3 * (k - 1))
        nc = number_of_cells(mesh)
        poro, perm = fill(0.2, nc), fill(5e-13, nc)
        flf = fine_fluid(p -> 0.8 * (2e3 / (max(p, 0.0) + 2e3))^2 + 0.2; cw = 4.3e-10)
        column(c) = cell_ijk(mesh, c)[1]
        # CO2 in columns 1–2 after 40 months, the strip of column 3 fine (:hyb)
        function updip(mode)
            isve = [mode == :ve || (mode == :hyb && column(c) != 3) for c in 1:nc]
            cols, _ = ve_columns(mesh; is_ve = isve)
            ve = ve_domain(mesh, cols; perm_x = perm, poro)
            fl = co2lab_ve_fluid(flf, mesh, cols; poro, perm_x = perm)
            # the injector in column 6: its VE block, or all its fine cells (a
            # simple well, its wellbore volume shrunk: the CO2 it starts with
            # would otherwise drain into the reservoir, 0.2% of the injected)
            wb = unique(cols.partition[c] for c in 1:nc if column(c) == 6)
            well = mode == :fine ?
                   setup_well(ve.domain, wb; name = :inj, WI = fill(1e-11, length(wb)),
                              simple_well = true, volume_multiplier = 1e-3) :
                   setup_ve_well(ve, only(wb); name = :inj)
            model = setup_co2lab_model(ve.domain; fluid = fl, wells = [well])
            on = setup_reservoir_forces(model; control = Dict(
                :inj => co2_injector_control(TotalRateTarget(q), fl)))
            off = setup_reservoir_forces(model;
                                         control = Dict(:inj => DisabledControl()))
            p0 = P .+ RHOW * g .* (ve.domain[:reference_depth] .- 1000.0)
            case = JutulCase(model, fill(30 * 86400.0, 40), [k <= 10 ? on : off
                                                             for k in 1:40];
                             state0 = setup_co2lab_state(model; pressure = p0,
                                                         forces = on),
                             parameters = setup_parameters(model))
            sim, cfg = setup_co2lab_simulator(case; info_level = -1)
            M = simulate_reservoir(case; simulator = sim,
                                   config = cfg).states[end][:TotalMasses][2, :]
            injected = 10 * 30 * 86400.0 * q * fl.rho_ref[2]
            @test abs(sum(M) - injected) / injected < 1e-4
            colb(b) = column(cols.cells[cols.pos[b]])
            return sum(M[b] for b in eachindex(M) if colb(b) < 3)
        end
        mve, mhyb, mfine = updip(:ve), updip(:hyb), updip(:fine)
        # 20 layers: 0.12 / 0.71 / 2.61 t
        @test mfine < mhyb < mve && mhyb > 0.2 * mve
    end

    @testset "exact checks P0, P2, P3" begin
        J = JutulCO2Lab
        g = Jutul.gravity_constant
        rw, rg = 0.2, 0.2
        curve = p -> 0.8 * (2e3 / (max(p, 0.0) + 2e3))^2 + 0.2
        # fringe, residual CO2 (hysteresis) and dissolution: every equation
        ff = FineScaleFluid(; rock_types = [RockCurves(; invPc3D = curve, kr3D = 2)],
                            residual = (rw, rg), default_pvt = nothing,
                            co2_rho_pvt = (4.3e-9, 1e7), wat_rho_pvt = (4.3e-10, 1e7),
                            dissolution = true, dis_max = 0.03)
        nt, dt = 24, 30 * 86400.0
        # fixed steps: the adaptive step control reacts to iteration counts,
        # so mirrored runs could take different steps (P3)
        kw = (info_level = -1, tol_cnv = 1e-10, tol_mb = 1e-12, timesteps = :none,
              max_timestep = Inf)
        # a closed run from CO2 placed in the blocks `co2` (no wells)
        function closed_run(mesh, cols, co2; kwarg...)
            nc = number_of_cells(mesh)
            poro, perm = fill(0.2, nc), fill(5e-13, nc)
            ve = ve_domain(mesh, cols; perm_x = perm, poro)
            fl = co2lab_ve_fluid(ff, mesh, cols; poro, perm_x = perm, dis_rate = 1e-10)
            model = setup_co2lab_model(ve.domain; fluid = fl)
            p0 = 1e7 .+ 1050.0 * g .* (ve.domain[:reference_depth] .- 1000.0)
            sg = [b in co2 ? 0.3 : 0.0 for b in 1:number_of_columns(cols)]
            case = JutulCase(model, fill(dt, nt), setup_reservoir_forces(model);
                             state0 = setup_co2lab_state(model; pressure = p0, sG = sg),
                             parameters = setup_parameters(model))
            return (model = model, ve = ve, case = case, kwarg...)
        end
        function simulate(r)
            sim, cfg = setup_co2lab_simulator(r.case; kw...)
            return simulate_reservoir(r.case; simulator = sim, config = cfg).states
        end

        # P0: on a VE grid of continuous faces the transition flux type gives
        # the stock flux exactly (its ordinary faces untouched)
        mesh = _cp_mesh_pillars(5, 8, (p, k) -> 1000.0 - 1.0 * p + 0.75 * (k - 1))
        cols, _ = ve_columns(mesh)
        r = closed_run(mesh, cols, (1, 2))
        @test J._ve_faces(r.ve.domain) === nothing
        st = simulate(r)
        N, zref = r.ve.domain[:neighbors], r.ve.domain[:reference_depth]
        faces = [J.VEFace(N[1, f], N[2, f], J._FACE_VE, NaN, (NaN, NaN, NaN, NaN),
                          (zref[N[1, f]], zref[N[2, f]]), (true, true))
                 for f in axes(N, 2)]
        J._set_transition_flux!(r.model[:Reservoir], faces)
        sf = simulate(r)
        @test maximum(st[end][:sG]) > 0.05 && maximum(st[end][:Rs]) > 1e-3
        for k in (:Pressure, :sG, :sGmax, :Rs)
            @test sf[end][k] == st[end][k]
        end

        # a stepped sloping grid (0.3 m per column) with a fault (offset 1.5
        # m, cells 0.6 m: the virtual cell) between columns 4 and 5 and a fine
        # pillar 2; CO2 placed in columns 1 and 3 migrates up-dip through the
        # fine pillar and towards the fault. `mirror` flips it left to right.
        nx = 6
        top(i, k) = 1000.0 - 0.3 * i + (i > 4 ? 1.5 : 0.0) + 0.6 * (k - 1)
        grid(mirror) = _cp_mesh(nx, 10, mirror ? (i, k) -> top(nx + 1 - i, k) : top)
        pillar(i, mirror) = mirror ? nx + 1 - i : i
        function hybrid(mirror)
            mesh = grid(mirror)
            cols, _ = ve_columns(mesh; is_ve = ve_except_pillars(mesh,
                                                [(pillar(2, mirror), 1)]))
            col(b) = cell_ijk(mesh, cols.cells[cols.pos[b]])[1]
            co2 = [b for b in 1:number_of_columns(cols) if cols.is_ve[b] &&
                   col(b) in pillar.((1, 3), mirror)]
            return closed_run(mesh, cols, co2; mesh = mesh, cols = cols)
        end

        # P2: brine and CO2 (free and dissolved) conserved across the
        # transitions, to the solver's tolerance
        a = hybrid(false)
        @test count(==(J._FACE_VE_DISCONTINUOUS), a.ve.domain[:face_category]) == 1
        @test count(==(J._FACE_VE_FINE), a.ve.domain[:face_category]) >= 20
        sa = simulate(a)
        M0 = sum(sa[1][:TotalMasses]; dims = 2)
        fine = .!a.cols.is_ve
        @test maximum(sa[end][:sG][fine]) > 0.01 && maximum(sa[end][:Rs][fine]) > 0.01
        for k in eachindex(sa)
            M = sum(sa[k][:TotalMasses]; dims = 2)
            @test maximum(abs.(M .- M0) ./ M0) < 1e-10
        end

        # P3: the mirrored grid gives the mirrored solution (each face's flux
        # independent of its orientation)
        b = hybrid(true)
        sb = simulate(b)
        # block of b holding the mirror image of a's block's top cell
        ijk_b = Dict(cell_ijk(b.mesh, c) => c for c in 1:number_of_cells(b.mesh))
        image = map(1:number_of_columns(a.cols)) do blk
            i, j, k = cell_ijk(a.mesh, a.cols.cells[a.cols.pos[blk]])
            b.cols.partition[ijk_b[(nx + 1 - i, j, k)]]
        end
        @test sort(image) == 1:number_of_columns(b.cols)
        for k in (:sG, :sGmax, :Rs)
            @test maximum(abs.(sb[end][k][image] .- sa[end][k])) < 1e-8
        end
        @test maximum(abs.(sb[end][:Pressure][image] .- sa[end][:Pressure])) < 0.1
    end

    @testset "fixed residual, capped (later item: co2lab's 'fixed residual')" begin
        J = JutulCO2Lab
        g = Jutul.gravity_constant
        mesh = _cp_mesh_pillars(4, 6, (p, k) -> 1000.0 - 2.0 * p + 1.0 * (k - 1))
        nc = number_of_cells(mesh)
        layer = [cell_ijk(mesh, c)[3] for c in 1:nc]
        perm = [1e-13 * (1 + 0.5 * sin(l)) for l in layer]
        poro = [0.15 + 0.01 * l for l in layer]
        rw, rg, pe = 0.2, 0.25, 2e3
        sand = p -> 0.8 * (pe / (max(p, 0.0) + pe))^2 + 0.2
        shale = p -> 0.8 * (4pe / max(p, 4pe))^2 + 0.2
        krg(s) = min(s, 1 - rw)^2
        krw(sw) = (max(sw - rw, 0.0) / (1 - rw))^2
        rt = [l in (3, 4) ? 2 : 1 for l in layer]
        cols, _ = ve_columns(mesh)
        dom = Dict(:column_height => fill(6.0, 4))
        drho, H = 300.0, 6.0
        fluid(law; kw = krw) = FineScaleFluid(;
            rock_types = [RockCurves(; invPc3D = sand, kr3D = krg, krw3D = kw),
                          RockCurves(; invPc3D = shale, kr3D = krg, krw3D = kw)],
            residual = (rw, rg), hysteresis = law, default_pvt = nothing)
        build(ff; types = rt) = (fl = co2lab_ve_fluid(ff, mesh, cols; poro, perm_x = perm,
                                                      rock_type = types);
                     J._vertical_model(fl.vertical, dom, fl))
        mf, me = build(fluid(FixedResidual())), build(fluid(EndpointScaling()))
        @test mf.law isa FixedResidual && mf.p_res[1] > 0 && mf.p_res[2] > 4pe

        # the depth rule by brute-force quadrature: per 1 m layer its rock type
        typ(z) = clamp(ceil(Int, z), 1, 6) in (3, 4) ? 2 : 1
        ph(z) = 0.15 + 0.01 * clamp(ceil(Int, z), 1, 6)
        kx(z) = 1e-13 * (1 + 0.5 * sin(clamp(ceil(Int, z), 1, 6)))
        sgz(z, hh) = hh > z ? 1 - (typ(z) == 1 ? sand : shale)(drho * g * (hh - z)) : 0.0
        loc(z, h, hm) = max(sgz(z, h), min(rg, sgz(z, hm)))
        mob(z, h, hm) = sgz(z, h) > rg ? sgz(z, h) : 0.0
        zs = range(0, H; length = 600_001)
        mid = zs[1:end-1] .+ step(zs) / 2
        quad(f) = sum(f(z) for z in mid) * step(zs)
        # (imbibed states hold at least the trapped CO2: 0.11 after 0.3)
        for (S, Smax) in ((0.15, 0.3), (0.12, 0.3), (0.3, 0.3), (0.03, 0.04))
            h, hm = J._column_heights(mf, 2, S, Smax, H, drho, rw, rg)
            @test h <= hm + 1e-12
            # the column's CO2 is S; relperms against the depth rule
            @test quad(z -> ph(z) * loc(z, h, hm)) / quad(ph) ≈ S rtol = 1e-5
            krW, krG = J._column_relperm(mf, 2, S, Smax, h, hm, H, drho, 0.0, 1.0, 1.0,
                                         rw, rg)
            K = quad(kx)
            @test quad(z -> kx(z) * krg(mob(z, h, hm))) / K ≈ krG rtol = 1e-4 atol = 1e-9
            @test quad(z -> kx(z) * krw(1 - loc(z, h, hm))) / K ≈ krW rtol = 1e-5
            # reconstruction: the cells' averages of the depth rule
            zb = 1.0:1.0:6.0
            s, _ = J._fine_saturations(mf, 2, collect(zb), h, hm, drho, rw, rg)
            @test s[2] ≈ quad(z -> 1.0 < z < 2.0 ? loc(z, h, hm) : 0.0) rtol = 1e-4
        end
        # drainage: as endpoint scaling
        @test all(J._column_heights(mf, 2, 0.2, 0.2, H, drho, rw, rg) .≈
                  J._column_heights(me, 2, 0.2, 0.2, H, drho, rw, rg))
        # the height grows with S at fixed S_max (monotone); below the trapped
        # CO2 (reached by dissolution only) linearly to 0
        hs = [J._column_heights(mf, 2, S, 0.3, H, drho, rw, rg)[1]
              for S in range(0.0, 0.3; length = 61)]
        @test issorted(hs) && hs[1] == 0

        # a curve whose fringe holds r_g from its tip on (a critical saturation
        # of r_g): co2lab's 'fixed residual', the
        # profile of h plus r_g over [h, h_max]
        tipped = p -> max(0.8 * (pe / (max(p, 0.0) + pe))^2 + 0.2 - rg, rw)
        mt = build(FineScaleFluid(; rock_types = [RockCurves(; invPc3D = tipped, kr3D = krg)],
                                  residual = (rw, rg), hysteresis = FixedResidual(),
                                  default_pvt = nothing); types = nothing)
        @test mt.p_res == [0.0]
        h, hm = J._column_heights(mt, 2, 0.12, 0.3, H, drho, rw, rg)
        I = J._interval_integrals(mt, 2, 0.0, H, h, hm, drho)
        pvh(a, b) = J._integral(mt.cells, 2, mt.cells.V, b) -
                    J._integral(mt.cells, 2, mt.cells.V, a)
        @test I.co2 ≈ I.free + rg * pvh(h, hm) rtol = 1e-10

        # below the trapped CO2 (S < the CO2 at h = 0): the derivative in S as
        # finite differences (the one in S_max omits the plateau end's, a
        # plain value: tiny here), and finite for a round-off trace of CO2
        # (subnormal S, S_max: the quotient rule's square underflowed to NaN)
        D = Jutul.ForwardDiff.Dual
        hS(S, Sm) = J._column_heights(mt, 2, S, Sm, H, drho, rw, rg)[1]
        ad = hS(D(0.05, 1.0, 0.0), D(0.3, 0.0, 1.0))
        fd(f, x, e) = (f(x + e) - f(x - e)) / 2e
        @test ad.partials[1] ≈ fd(S -> hS(S, 0.3), 0.05, 1e-6) rtol = 1e-5
        tiny = hS(D(5e-320, 1.0, 0.0), D(6e-320, 0.0, 1.0))
        @test all(isfinite, tiny.partials) && tiny.value < 1e-300

        # the trapping inventory's split: immobile the residual min(r_g, s_max)
        # (what stays when the mobile CO2 has gone), mobile the rest; the
        # parts sum to S; a trap 2 m deep and a 0.5 m relief in column 2
        Sv, Smv = [0.15, 0.12, 0.3, 0.03], [0.3, 0.3, 0.3, 0.04]
        hh = [J._column_heights(mf, i, Sv[i], Smv[i], H, drho, rw, rg) for i in 1:4]
        col = (H = fill(H, 4), dh = [0.0, 0.5, 0.0, 0.0], zt = [0.0, 2.0, 1.0, 0.0])
        c = J._category_saturations(mf, Sv, Smv, first.(hh), last.(hh), col,
                                    fill(drho, 4), rw, rg)
        parts = c.immob_sub .+ c.immob_struct .+ c.immob_free .+ c.mobile_sub .+
                c.mobile_struct .+ c.mobile_free
        @test parts ≈ Sv rtol = 1e-8
        @test all(c.res .<= c.immob_free .+ 1e-12)
        h2, hm2 = hh[2]
        res2(z) = min(rg, sgz(z, hm2))
        @test c.immob_struct[2] ≈ quad(z -> 0.5 < z < 2.0 ? ph(z) * res2(z) : 0.0) /
                                 quad(ph) rtol = 1e-4
        @test c.mobile_struct[2] ≈
              quad(z -> 0.5 < z < 2.0 ? ph(z) * (loc(z, h2, hm2) - res2(z)) : 0.0) /
              quad(ph) rtol = 1e-4
        # co2lab's split (endpoint scaling's trapped amount) on request
        cc = J._category_saturations(mf, Sv, Smv, first.(hh), last.(hh), col,
                                     fill(drho, 4), rw, rg; co2lab_compatible = true)
        @test cc.immob_free .+ cc.immob_struct .+ cc.immob_sub ≈ Smv .* rg ./ (1 - rw)

        # derivatives (one Newton step with AD) against central differences
        f(x) = (hh = J._column_heights(mf, 2, x[1], x[2], H, x[3], rw, rg);
                [hh..., J._column_relperm(mf, 2, x[1], x[2], hh..., H, x[3], 0.0, 1.0,
                                          1.0, rw, rg)...])
        for x0 in ([0.15, 0.3, 300.0], [0.25, 0.3, 280.0])
            jac = Jutul.ForwardDiff.jacobian(f, x0)
            e = [1e-7, 1e-7, 1e-4] .* x0
            fd = reduce(hcat, [(f(x0 .+ e[j] .* (1:3 .== j)) .-
                                f(x0 .- e[j] .* (1:3 .== j))) ./ (2e[j]) for j in 1:3])
            @test jac ≈ fd rtol = 1e-4
        end

        # runs, VE and hybrid (a fine pillar), injection then migration: CO2
        # conserved, imbibition reached, and more CO2 left behind trapped than
        # with endpoint scaling (r_g where the plume swept, not r_g·s_max/(1 − r_w))
        slope = _cp_mesh_pillars(6, 6, (p, k) -> 1000.0 - 2.0 * p + 1.0 * (k - 1))
        ns = number_of_cells(slope)
        ps, ks = fill(0.2, ns), fill(5e-13, ns)
        function run(law, is_ve)
            c6, _ = ve_columns(slope; is_ve)
            ve = ve_domain(slope, c6; perm_x = ks, poro = ps)
            ff = FineScaleFluid(; rock_types = [RockCurves(; invPc3D = sand, kr3D = krg)],
                                residual = (rw, rg), hysteresis = law,
                                default_pvt = nothing, wat_rho_pvt = (4.3e-10, 1e7),
                                co2_rho_pvt = (4.3e-9, 1e7))
            fl = co2lab_ve_fluid(ff, slope, c6; poro = ps, perm_x = ks)
            wb = c6.partition[findfirst(c -> cell_ijk(slope, c) == (1, 1, 1), 1:ns)]
            well = c6.is_ve[wb] ? setup_ve_well(ve, wb; name = :inj) :
                   setup_well(ve.domain, [wb]; name = :inj, WI = 1e-11)
            model = setup_co2lab_model(ve.domain; fluid = fl, wells = [well])
            q = 2e-4
            on = setup_reservoir_forces(model; control = Dict(
                :inj => co2_injector_control(TotalRateTarget(q), fl)))
            off = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()))
            p0 = 1e7 .+ 1050.0 * g .* (ve.domain[:reference_depth] .- 990.0)
            nt = 24
            case = JutulCase(model, fill(30 * 86400.0, nt), [k <= 6 ? on : off
                                                             for k in 1:nt];
                             state0 = setup_co2lab_state(model; pressure = p0,
                                                         forces = on),
                             parameters = setup_parameters(model))
            sim, cfg = setup_co2lab_simulator(case; info_level = -1)
            st = simulate_reservoir(case; simulator = sim, config = cfg).states[end]
            injected = 6 * 30 * 86400.0 * q * fl.rho_ref[2]
            @test abs(sum(st[:TotalMasses][2, :]) - injected) / injected < 1e-4
            return st
        end
        allve = fill(true, ns)
        hyb = ve_except_pillars(slope, [(3, 1)])
        sf, se = run(FixedResidual(), allve), run(EndpointScaling(), allve)
        @test any(sf[:sG] .< sf[:sGmax] .- 0.01)
        @test sum(sf[:sG] .* (sf[:sG] .< sf[:sGmax])) > sum(se[:sG] .* (se[:sG] .< se[:sGmax]))
        sh = run(FixedResidual(), hyb)
        @test any(sh[:sG] .< sh[:sGmax] .- 0.01)

        # fine cells: the drainage curve until trapped at min(r_g, s_max)
        @test J._fine_kr_saturation(FixedResidual(), 0.3, 0.5, rw, rg) == 0.3
        @test J._fine_kr_saturation(FixedResidual(), 0.25, 0.5, rw, rg) == 0.0
        @test J._fine_kr_saturation(FixedResidual(), 0.1, 0.1, rw, rg) == 0.1
        @test J._fine_pc_saturation(FixedResidual(), 0.1, 0.5, rw, rg) == 0.1
    end

    @testset "Newton updates stopped at relperm kinks in fine cells" begin
        J = JutulCO2Lab
        mesh = _cp_mesh(3, 4, (i, k) -> 1000.0 + 1.0 * (k - 1))
        nc = number_of_cells(mesh)
        perm, poro = fill(1e-13, nc), fill(0.2, nc)
        # a critical saturation of 0.15 (= r_g, one kink), or of 0.1 (two)
        crit(sc) = s -> s <= sc ? 0.0 : (s - sc) / (0.8 - sc)
        ff(sc) = FineScaleFluid(; rock_types = [RockCurves(;
                                    invPc3D = p -> max(0.85 - p / 1e5, 0.2),
                                    kr3D = crit(sc))],
                                residual = (0.2, 0.15), hysteresis = FixedResidual(),
                                default_pvt = nothing)
        @test J._saturation_kinks(ff(0.15)) ≈ [[0.15]]
        @test J._saturation_kinks(ff(0.1)) ≈ [[0.1, 0.15]]
        @test J._crossed_kink([0.1, 0.15], 0.05, 0.2) == 0.1
        @test J._crossed_kink([0.1, 0.15], 0.2, 0.05) == 0.15
        @test isnothing(J._crossed_kink([0.15], 0.15, 0.3))

        # fine cells get the safeguard, an all-VE domain keeps the plain unknown
        sG(cols) = (fl = co2lab_ve_fluid(ff(0.1), mesh, cols; poro, perm_x = perm);
                    ve = ve_domain(mesh, cols; perm_x = perm, poro);
                    setup_co2lab_model(ve.domain; fluid = fl)[:Reservoir])
        fine = sG(first(ve_columns(mesh; is_ve = fill(false, nc))))
        @test sG(first(ve_columns(mesh))).primary_variables[:sG] isa J.GasSaturation
        v = fine.primary_variables[:sG]
        @test v isa J.FineGasSaturation && all(==(1), v.rock_type)

        # an update across a kink stops at the nearest one; others as Jutul's
        st = Dict{Symbol, Any}(:sG => [0.05; 0.12; 0.3; 0.2; fill(0.5, nc - 4)])
        dx = [0.15; 0.1; -0.05; -0.18; zeros(nc - 4)]
        Jutul.update_primary_variable!(st, v, :sG, fine, dx, 1.0)
        @test st[:sG][1:4] ≈ [0.1, 0.15, 0.25, 0.15] && all(==(0.5), st[:sG][5:end])
    end

    @testset "endpoint-scaling hysteresis in fine 3D cells" begin
        J = JutulCO2Lab
        rw, rg = 0.2, 0.25
        krg(s) = min(s, 1 - rw)^2
        imb(s) = krg((1 - rw) * max(s - rg, 0.0) / (1 - rw - rg))
        s = collect(range(0.0, 1.0; length = 4001))
        tab(f) = JutulDarcy.PhaseRelativePermeability(s, f.(s); label = :g)
        H = EndpointScalingHysteresis(rw, rg)
        # the relperm at co2lab's free saturation, below and at the largest
        # maximum (JutulDarcy takes the imbibition curve there)
        for (sg, smax) in ((0.3, 0.5), (0.1, 0.4), (0.5, 0.8), (0.3, 0.3))
            kr = JutulDarcy.kr_hysteresis(H, tab(krg), tab(imb), sg, smax)
            @test kr ≈ krg(J._free_saturation(sg, smax, rw, rg)) atol = 1e-6
        end

        # a column of 10 cells: CO2 injected at the bottom, then brine; cells
        # with only residual CO2 have the capillary pressure of zero free CO2
        mesh = reservoir_mesh((1, 1, 10), (10.0, 10.0, 10.0))
        domain = reservoir_domain(mesh; permeability = 1e-13, porosity = 0.2)
        sys = ImmiscibleSystem((AqueousPhase(), VaporPhase());
                               reference_densities = (1000.0, 700.0))
        well = setup_well(domain, [10]; name = :inj)
        model = setup_reservoir_model(domain, sys; wells = [well])
        pc(sg) = 2e3 * sg
        kr = JutulDarcy.ReservoirRelativePermeabilities(
            w = JutulDarcy.PhaseRelativePermeability(s, (x -> max(x - rw, 0.0) /
                                                     (1 - rw)).(s); label = :w),
            g = (tab(krg), tab(imb)), hysteresis_g = H)
        replace_variables!(model; RelativePermeabilities = kr)
        JutulDarcy.add_relperm_parameters!(model)
        set_secondary_variables!(model[:Reservoir];
                                 CapillaryPressure = FreeSaturationCapillaryPressure(pc,
                                                                                    rw, rg))
        push!(model[:Reservoir].output_variables, :MaxSaturations, :CapillaryPressure)
        z = domain[:cell_centroids][3, :]
        state0 = setup_reservoir_state(model; Pressure = 1e7 .+ 1e4 .* z,
                                       Saturations = [1.0, 0.0])
        ctrl(f) = Dict(:inj => InjectorControl(TotalRateTarget(2e-5), f;
                                               density = f[1] > 0 ? 1000.0 : 700.0))
        forces = [setup_reservoir_forces(model; control = ctrl(f))
                  for f in ([0.0, 1.0], [1.0, 0.0])]
        dt = fill(10 * 86400.0, 20)
        res = simulate_reservoir(state0, model, dt; forces = [k <= 10 ? forces[1] :
                                                             forces[2] for k in 1:20],
                                 info_level = -1)
        st = res.states[end]
        sg, smax = st[:Saturations][2, :], st[:MaxSaturations][2, :]
        free = J._free_saturation.(sg, smax, rw, rg)
        @test st[:CapillaryPressure][1, :] ≈ pc.(free)
        # residual CO2 left behind by the brine, holding no capillary pressure
        trapped = findall(i -> sg[i] > 1e-3 && free[i] < 1e-8, eachindex(sg))
        @test !isempty(trapped)
        @test all(abs.(st[:CapillaryPressure][1, trapped]) .< 1e-3)
    end

    @testset "tabulated PVT" begin
        J = JutulCO2Lab
        # JutulDarcy's CSP11 CO2 viscosity is NIST-consistent, as co2lab's
        # (co2lab's values evaluated under Octave)
        dump = joinpath(refdir, "co2lab_co2_tables.mat")
        if isfile(dump)
            d = matread(dump)
            i, j = findfirst(==(20e6), vec(d["p"])), findfirst(==(343.0), vec(d["T"]))
            jd = jutuldarcy_property_tables()
            @test abs(jd.viscosity[2](20e6, 343.0) / d["mu"][i, j] - 1) < 0.02
            @test abs(jd.density[2](20e6, 343.0) / d["rho"][i, j] - 1) < 1e-3
        else
            @warn "co2lab's table dump not found"
        end

        # the fluid's defaults: JutulDarcy's tables, surface densities, and
        # co2lab's default dissolution capacity as a mass ratio
        fl = co2lab_ve_fluid()
        @test needs_temperature(fl) && all(x -> x isa PhaseTable, (fl.rho_pvt...,
                                                                   fl.mu_pvt...))
        @test fl.rho_ref == (998.207150, 1.868048)
        @test fl.dis_max * fl.rho_ref[2] / fl.rho_ref[1] ≈ 0.03 * 760 / 1100
        @test co2lab_ve_fluid(; default_pvt = nothing).dis_max == 0.03
        @test co2lab_ve_fluid(; dis_max_mass = 0.05).dis_max ≈
              0.05 * 998.20715 / 1.868048
        # reference densities at reservoir conditions, from the tables
        rc = co2lab_ve_fluid(; reference_conditions = (20e6, 343.15))
        @test rc.rho_ref[2] ≈ rc.rho_pvt[2](20e6, 343.15)
        # a reference value without its property model is refused
        @test_throws ErrorException co2lab_ve_fluid(; co2_rho_ref = 700.0)
        @test_throws ErrorException co2lab_ve_fluid(; dis_max = 0.03,
                                                    dis_max_mass = 0.02)
        # salinity needs the Salo correlations
        @test_throws ErrorException jutuldarcy_property_tables(;
            salt_names = ["NaCl"], salt_mole_fractions = [0.01])
        # mixed: tabulated CO2, co2lab's constant brine
        mix = co2lab_ve_fluid(; wat_rho_pvt = nothing, wat_mu_pvt = nothing)
        @test mix.rho_ref[1] == 1100.0 && mix.rho_pvt[1] == (0.0, 0.0)

        # a run with a geothermal gradient: temperatures per column and in the
        # well, CO2 conserved, and the CO2 density following the temperature
        mesh = _cp_mesh(4, 3, (i, k) -> 1500.0 + 4.0 * i + 5.0 * (k - 1))
        cols, _ = ve_columns(mesh)
        grad = geothermal_gradient(T_ref = 283.15, gradient = 0.03)
        ve = ve_domain(mesh, cols; perm_x = fill(1e-13, 12), poro = fill(0.2, 12),
                       temperature = grad)
        @test ve.domain[:temperature] ≈ grad.(ve.geometry.z_top)
        well = setup_ve_well(ve, 4; name = :inj)
        model = setup_co2lab_model(ve.domain; wells = [well], fluid = fl)
        prm = setup_parameters(model)
        @test prm[:inj][:Temperature] ≈ fill(ve.domain[:temperature][4], 1)
        q = 0.5 / fl.rho_ref[2]
        on = setup_reservoir_forces(model; control = Dict(:inj => co2_injector_control(
                                        TotalRateTarget(q), fl)))
        p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = 1000.0)
        case = JutulCase(model, fill(10 * 86400.0, 3), fill(on, 3); parameters = prm,
                         state0 = setup_co2lab_state(model; pressure = p0, forces = on))
        sim, cfg = setup_co2lab_simulator(case; info_level = -1)
        res = simulate_reservoir(case; simulator = sim, config = cfg)
        injected = 3 * 10 * 86400.0 * 0.5
        # the CO2 not in the reservoir is in the wellbore, which fills further
        # as the pressure rises
        missing_co2 = injected - sum(res.states[end][:TotalMasses][2, :])
        @test 0 <= missing_co2 <= sum(res.result.states[end][:inj][:TotalMasses][2, :])
        b(T) = J._shrinkage(fl, 20e6, T)[2]
        @test b(320.0) > b(360.0)
    end

    @testset "inventory split of a fringe column" begin
        J = JutulCO2Lab
        fl = co2lab_ve_fluid(; default_pvt = nothing,
                             vertical = CapillaryFringeTable(; invPc3D = (0.4, 0.5),
                                                               kr3D = 2),
                             residual = (0.2, 0.2), co2_rho_ref = 700.0,
                             wat_rho_ref = 1020.0)
        H, drho = 20.0, 320.0
        m = J._vertical_model(fl.vertical, Dict(:column_height => [H]), fl)
        S, Smax = [0.3], [0.45]
        h, hmax = J._column_heights(m, 1, S[1], Smax[1], H, drho, 0.2, 0.2)
        col = (H = [H], zt = [5.0], dh = [1.0])
        split(c) = J._category_saturations(m, S, Smax, [h], [hmax], col, [drho], 0.2,
                                           0.2; co2lab_compatible = c)
        p, c = split(false), split(true)
        # co2lab's switch is the sharp interface's split by height
        @test c == J._height_split(S, Smax, [h], [hmax], col, 0.2, 0.2)
        # the parts add up, and the profile puts more CO2 under the caprock
        @test p.immob_sub .+ p.immob_struct .+ p.immob_free ≈ Smax .* 0.2 ./ 0.8
        @test p.mobile_sub .+ p.mobile_struct .+ p.mobile_free ≈ S .- Smax .* 0.25
        @test all(p.mobile_struct .> c.mobile_struct)
        @test all(p.immob_struct .> c.immob_struct)
        @test all(0 .<= p.res .<= p.immob_free)
    end

    @testset "CO2 Storage Atlas maps vs co2lab" begin
        f = joinpath(refdir, "atlas_maps.mat")
        if !isfile(f) || !isdir(atlas_data_directory())
            @info "Skipping: Atlas data or MRST reference dump $f missing"
        else
            rep = _atlas_vs_mrst(matread(f))
            passed(rep) || display(rep)
            @test passed(rep)
        end

        # co2lab's row handling: its map is ours one row further north, with
        # the southernmost row empty
        if isdir(atlas_data_directory())
            path = first(atlas_map_files()).path
            a = read_aai_grid(path).data
            b = read_aai_grid(path; co2lab_rows = true).data
            @test all(isnan, b[:, 1]) && isequal(b[:, 2:end], a[:, 1:end-1])
            @test all(isnan, a[:, end])       # the Atlas files' first row is empty
        end

        # names and variants as co2lab builds them; average rock by name
        files = atlas_map_files()
        @test any(f -> f.name == "Huginfmeast" && f.variant == :thickness, files)
        @test any(f -> f.name == "Jurassicmid" && f.variant == :top, files)
        @test atlas_average_rock("Tiljefm").ntg == 0.30
        @test isnan(atlas_average_rock("Rorfm").perm)
        @test length(atlas_formation_names()) == 31
        @test_throws ArgumentError atlas_formation_names(:arctic)
    end

    @testset "CO2 Storage Atlas grid sections" begin
        f = joinpath(refdir, "atlas_grids.mat")
        if !isfile(f) || !isdir(atlas_data_directory())
            @info "Skipping: Atlas data or MRST reference dump $f missing"
        else
            rep = _atlas_grids_vs_mrst(matread(f))
            passed(rep) || display(rep)
            @test passed(rep)
        end

        # bilinear interpolation reproduces a bilinear function
        xs, ys = [0.0, 1.0, 3.0], [10.0, 12.0]
        fb(x, y) = 1 + 2x - y + 0.5x * y
        Z = [fb(x, y) for y in ys, x in xs]
        @test JutulCO2Lab._interp2(xs, ys, Z, 2.2, 11.3) ≈ fb(2.2, 11.3)
        @test isnan(JutulCO2Lab._interp2(xs, ys, Z, 3.1, 11.0))
        @test JutulCO2Lab._linspace(0.0, 1.0, 5) == [0.0, 0.25, 0.5, 0.75, 1.0]

        if isdir(atlas_data_directory())
            a = atlas_grid_section("Johansenfm"; coarsening = 3, co2lab_compatible = true)
            b = atlas_grid_section("Johansenfm"; coarsening = 3)
            # values at their cell centres: half the file's cell size (200 m)
            ca, cb = reshape(a.grid["COORD"], 6, :), reshape(b.grid["COORD"], 6, :)
            @test cb[1:2, 1] ≈ ca[1:2, 1] .+ 100.0
            # every active cell has four defined pillars
            nx, ny, _ = b.grid["cartDims"]
            ztop = reshape(cb[3, :], nx + 1, ny + 1)
            @test all(b.grid["ACTNUM"][i, j, 1] <= all(ztop[i:i+1, j:j+1] .> 0)
                      for i in 1:nx, j in 1:ny)
            m = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(b.grid)
            @test number_of_cells(m) == count(b.grid["ACTNUM"])

            # Jurassic formations use the Jurassic top; rock maps on all layers
            @test count(atlas_grid_section("Huginfmwest"; coarsening = 4).grid["ACTNUM"]) > 0
            s = atlas_grid_section("Stofm"; coarsening = 4, nz = 2).grid
            @test size(s["PORO"]) == s["cartDims"]
            @test isequal(s["PORO"][:, :, 1], s["PORO"][:, :, 2])
            @test_throws ArgumentError atlas_grid_section("Paleocene")
        end
    end

    @testset "CO2 Storage Atlas VE domains" begin
        f = joinpath(refdir, "atlas_formations.mat")
        if !isfile(f) || !isdir(atlas_data_directory())
            @info "Skipping: Atlas data or MRST reference dump $f missing"
        else
            rep = _atlas_domains_vs_mrst(matread(f))
            passed(rep) || display(rep)
            @test passed(rep)
        end

        if isdir(atlas_data_directory())
            # one connected piece; with rock maps, no NaN left
            o = atlas_ve_domain("Stofm"; coarsening = 4, layers = 2)
            N = o.ve.domain[:neighbors]
            g = JutulCO2Lab.Graphs.SimpleGraph(number_of_columns(o.columns))
            foreach(f -> JutulCO2Lab.Graphs.add_edge!(g, N[1, f], N[2, f]), axes(N, 2))
            @test JutulCO2Lab.Graphs.is_connected(g)
            @test !any(isnan, o.ve.rock.perm) && !any(isnan, o.rock.ntg)
            @test all(==(2), diff(o.columns.pos))
        end
    end

    @testset "All CO2 Storage Atlas formations" begin
        f = joinpath(refdir, "atlas_all.mat")
        if !isfile(f) || !isdir(atlas_data_directory())
            @info "Skipping: Atlas data or MRST reference dump $f missing"
        else
            rep = _atlas_all_vs_mrst(matread(f))
            passed(rep) || display(rep)
            @test passed(rep)
        end
        if isdir(atlas_data_directory())
            @test "Johansenfm" in atlas_grid_formations(:north_sea)
            @test !("Paleocene" in atlas_grid_formations())      # top map only
            @test_throws ArgumentError atlas_ve_domain("Brentgrp"; coarsening = 10)
        end
    end

    @testset "Sleipner benchmark" begin
        # without the data: an error naming the directory, not a failure later
        withenv("CO2LAB_ATLAS_DATA" => mktempdir()) do
            err = try sleipner_grid_section(); nothing catch e; e end
            @test err isa ArgumentError && occursin("no longer publicly", err.msg)
        end
        # with the data: the grid section only (the 43-layer VE model is
        # compared outside the test suite)
        if isdir(sleipner_data_directory())
            g = sleipner_grid_section()
            @test g["cartDims"] == (65, 119, 43)
            @test count(g["ACTNUM"]) == 262990          # co2lab's G.cells.num
            c = reshape(g["COORD"], 3, :)
            @test extrema(c[1, :]) == (436914.0, 440114.0)   # MAPAXES: UTM zone 31
        end
    end

    @testset "top surface and contour lines" begin
        # sloping aquifer, 5 layers: one quadrilateral per column
        f = joinpath(refdir, "sloping_aquifer.mat")
        if isfile(f)
            ve, cols, _ = _ve_case(matread(f))
            m = ve_mesh_from_mrst(matread(f)["G3"])
            s = ve_top_surface(m, cols)
            @test length(s.polygons) == number_of_columns(cols)
            @test all(p -> length(p) == 4, s.polygons)
            # a 100 x 1 strip: outline = 2 long sides of 100 edges + 2 ends
            @test length(s.outline) == 202

            # a field linear in x is contoured exactly at x = level
            xa, xb = extrema(s.x)
            levels = [xa + 0.3017 * (xb - xa), xa + 0.7123 * (xb - xa)]
            lines = ve_contour_lines(s, s.x, levels)
            @test all(!isempty, lines)
            @test all(abs(seg[1] - levels[1]) < 1e-8 * xb &&
                      abs(seg[3] - levels[1]) < 1e-8 * xb for seg in lines[1])
            # per-column values; NaN columns are skipped
            v = collect(1.0:number_of_columns(cols)); v[1:10] .= NaN
            @test !isempty(only(ve_contour_lines(s, v, [50.5])))
            @test_throws ArgumentError ve_contour_lines(s, [1.0, 2.0], [1.0])
        end

        # an Atlas map as a lattice surface, NaN outside the formation
        if isdir(atlas_data_directory())
            mp = only(filter(m -> m.variant == :top,
                             read_atlas_maps(["Johansenfm"]; coarsening = 4)))
            ms = JutulCO2Lab._atlas_map_surface(mp)
            @test length(ms.polygons) == prod(size(mp.data) .- 1)
            seg = only(ve_contour_lines(ms, ms.z, [2500.0]))
            x0, x1 = extrema(ms.x)
            @test !isempty(seg) && all(x0 <= t[1] <= x1 for t in seg)
            @test length(atlas_coastline()) > 100 && length(atlas_wells().x) > 1000

            # the coastline, converted from UTM 33 to the Atlas's UTM 31, passes
            # within 1 km of Lindesnes (57.9826 N, 7.0475 E)
            x, y = JutulCO2Lab._utm_forward(57.9826, 7.0475, 31)
            c = reduce(vcat, atlas_coastline())
            @test minimum(filter(!isnan, hypot.(c[:, 1] .- x, c[:, 2] .- y))) < 1000
        end

        # UTM: round trip to below 1 mm, zone change consistent
        for (lat, lon) in ((57.98, 7.05), (71.17, 25.78), (60.4, -1.0)), z in (31, 33)
            la, lo = JutulCO2Lab._utm_inverse(JutulCO2Lab._utm_forward(lat, lon, z)..., z)
            @test hypot(la - lat, (lo - lon) * cosd(lat)) * 111e3 < 1e-3
        end
        let (x33, y33) = JutulCO2Lab._utm_forward(60.4, 5.3, 33)
            @test all(isapprox.(JutulCO2Lab._utm_rezone(x33, y33, 33, 31),
                                JutulCO2Lab._utm_forward(60.4, 5.3, 31); atol = 1e-3))
        end
    end

    @testset "sheets: trap analysis per internal layer" begin
        # a sloping aquifer, 12 x 6 x 6, with a sealed barrier between layers 3
        # and 4 under pillars 4-8 x 2-5 (the first trap's dome)
        grid = sloping_aquifer_grid(; nx = 12, ny = 6, nz = 6, length = 1200.0,
                                    width = 600.0, thickness = 12.0, slope = 0.01,
                                    traps = ((0.5, 0.4, 6.0),))
        m = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
        ijk = [cell_ijk(m, c) for c in 1:number_of_cells(m)]
        patch(i, j) = 4 <= i <= 8 && 2 <= j <= 5
        N = get_neighborship(m)
        layer(f, k) = (a = ijk[N[1, f]]; b = ijk[N[2, f]];
                       a[1:2] == b[1:2] && patch(a[1], a[2]) && minmax(a[3], b[3]) == k)
        barrier = [f for f in axes(N, 2) if layer(f, (3, 4))]
        plain, _ = ve_columns(m)
        cols, _ = ve_columns(m; transmissibility_multiplier = _multipliers(m, barrier),
                             below_holes = :stack)
        sh = ve_sheets(m, cols)
        @test sh.nsheets == 2 && sh.top == [true, false]
        @test count(==(2), sh.sheet) == 20 && count(==(1), sh.sheet) == 72
        @test ve_sheets(m, plain).nsheets == 1

        # the top sheet: the plain grid's surface (only the column heights differ)
        s0, s1 = trap_surface(m, plain), trap_surface(m, cols)
        @test s0.blocks == 1:72 && s1.blocks == findall(==(1), sh.sheet)
        @test all(getproperty(s0, k) == getproperty(s1, k)
                  for k in (:x, :y, :z, :polygons, :edges, :edge_columns, :sealed, :zc))
        r0, r1 = trap_analysis(s0), trap_analysis(s1)
        @test r0.traps == r1.traps && r0.trap_z == r1.trap_z

        # the sheet under the barrier: its own trap, every rim edge open, CO2
        # leaving it rises into the neighbouring pillar's column (the top sheet)
        s2 = trap_surface(m, cols; sheet = 2)
        @test length(s2.blocks) == 20 && !any(s2.sealed)
        @test all(cols.column_ij[b][1] in 4:8 for b in s2.blocks)
        @test length(trap_analysis(s2).trap_z) == 1
        t = sheet_spill_targets(m, cols, s2)
        @test length(t) == 2 * (5 + 4) && all(x -> x.target_sheet == 1, t)
        adjacent(a, b) = sum(abs, cols.column_ij[a] .- cols.column_ij[b]) == 1
        @test all(x -> adjacent(x.block, x.target) && sh.sheet[x.block] == 2, t)
        @test_throws ArgumentError trap_surface(m, cols; sheet = 3)

        # the barrier as a fine layer (layers 3-4 fine) and a fine pillar beside
        # the patch: the same lower sheet; CO2 over that rim edge stays fine
        isve = [!(patch(x[1], x[2]) && x[3] in (3, 4)) && x[1:2] != (9, 3) for x in ijk]
        hc, _ = ve_columns(m; is_ve = isve)
        hs = ve_sheets(m, hc)
        @test hs.nsheets == 2 && count(==(2), hs.sheet) == 20 &&
              count(==(0), hs.sheet) == 46
        ht = sheet_spill_targets(m, hc, trap_surface(m, hc; sheet = 2))
        fine = only(filter(x -> hc.column_ij[x.block] == (8, 3) &&
                                hc.column_ij[x.target] == (9, 3), ht))
        @test fine.target_sheet == 0 && !hc.is_ve[fine.target] &&
              hc.cells[hc.pos[fine.target]] == findfirst(==((9, 3, 1)), ijk)

        # a barrier across the whole width: the columns under it touch the
        # outer boundary; their sides' edges lie at their tops (the barrier)
        across = [f for f in axes(N, 2) if (a = ijk[N[1, f]]; b = ijk[N[2, f]];
                  a[1:2] == b[1:2] && 4 <= a[1] <= 8 && minmax(a[3], b[3]) == (3, 4))]
        fc, _ = ve_columns(m; transmissibility_multiplier = _multipliers(m, across),
                           below_holes = :stack)
        nc = number_of_cells(m)
        fv = ve_domain(m, fc; perm_x = fill(1e-13, nc), poro = fill(0.2, nc))
        fs, _ = ve_boundary_sides(m, fc, fv.domain)
        lower = findall(.!JutulCO2Lab._first_blocks(fc))
        # those in the outer rows (j = 1, 6) have sides on the boundary
        outer = filter(b -> fc.column_ij[b][2] in (1, 6), lower)
        @test length(lower) == 30 && length(outer) == 10 && all(in(fs.column), outer)
        # each such side's edge 6 m (three cells) below the same side of the
        # column above
        edge(b, d) = only(fs.edge_depth[i] for i in eachindex(fs.column)
                          if fs.column[i] == b && fs.direction[i] == d)
        side_of(b) = fc.column_ij[b][2] == 1 ? 3 : 4
        @test all(isapprox(edge(b, side_of(b)), edge(b - 1, side_of(b)) + 6.0;
                           atol = 1e-9) for b in outer)

        # trap heights for the inventory: per sheet, mapped to the columns
        H = ve_column_geometry(m, cols).H
        zt = JutulCO2Lab._trap_heights(H, [s1, s2], [r1, trap_analysis(s2)])
        @test any(>(0), zt[s2.blocks]) && any(>(0), zt[s1.blocks])
        @test zt == JutulCO2Lab._trap_heights(H, [s2, s1], [trap_analysis(s2), r1])
    end

    @testset "trapping analysis vs co2lab-spillpoint" begin
        f = joinpath(refdir, "traps.mat")
        if !isfile(f) || !isdir(atlas_data_directory())
            @info "Skipping: Atlas data or MRST reference dump $f missing"
        else
            d = matread(f)
            for rep in (_trap_surface_vs_mrst(d), _traps_vs_mrst(d, :cell),
                        _trap_quantities_vs_mrst(d, :cell), _traps_vs_mrst(d, :node),
                        _trap_quantities_vs_mrst(d, :node), _closed_edges_vs_mrst(d))
                passed(rep) || display(rep)
                @test passed(rep)
            end
        end
        @test_throws ArgumentError trap_analysis((;); method = :other)
        @test_throws ArgumentError trap_analysis((;); closed_fault_edges = [(1, 2)])

        # sealing faults (default) close the surface; open ones leave a complete
        # offset an open boundary, as co2lab: the west block's fault trap
        mesh, cols = _faulted_mesh(40.0)
        sealing = trap_surface(mesh, cols)
        stitched = trap_surface(mesh, cols; faults = :open)
        @test count(sealing.sealed) == 14 && !any(stitched.sealed)
        @test all(==(0), sealing.edge_columns[2, sealing.sealed])
        crest = 5 + 3 * 10
        for method in (:cell, :node)
            @test trap_analysis(sealing; method).traps[crest] == 1
            @test maximum(trap_analysis(stitched; method).traps; init = 0) == 0
        end
        @test_throws ArgumentError trap_surface(mesh, cols; faults = :leaky)

        # IGEMS surfaces: the file's orientation by default, co2lab's transposed
        if isdir(atlas_data_directory())
            a = trap_surface(igems_mesh("OSS", 1; coarsening = 4)...)
            b = trap_surface(igems_mesh("OSS", 1; coarsening = 4,
                                        co2lab_compatible = true)...)
            @test Set(zip(a.x, a.y, a.z)) == Set(zip(b.y, b.x, b.z))
            @test_throws ArgumentError igems_surface_path("OSS", 1000)
        end
    end

    @testset "capacity estimate" begin
        f = joinpath(refdir, "capacity.mat")
        if isfile(f)
            rep = _capacity_vs_mrst(matread(f))
            passed(rep) || display(rep)
            @test passed(rep)
        else
            @info "Skipping: MRST reference dump $f missing"
        end

        # the synthetic caprock, a constant CO2 density: the parts add up to the
        # pore volume, the trap parts to the trap volumes
        s = trap_surface(_sine_mesh()...)
        res = trap_analysis(s; method = :node)
        n, poro, rho = length(s.zc), 0.25, 700.0
        p, T = fill(10e6, n), fill(320.0, n)
        dense(p, T) = rho
        a = trapping_capacity(s, res; p, T, res_gas = 1.0, res_water = 0.0,
                              dis_max = 0.0, poro, co2_density = dense)
        @test sum(a.structural) ≈ rho * sum(trap_volumes(s, res; poro = fill(poro, n)))
        @test sum(a.structural) + sum(a.residual) ≈ rho * poro * sum(s.area .* s.H)
        @test all(iszero, a.dissolved)
        b = trapping_capacity(s, res; p, T, res_gas = 0.2, res_water = 0.1,
                              dis_max = 50.0, poro, co2_density = dense)
        pv_trap = s.area .* b.trap_height .* poro
        @test b.structural ≈ pv_trap .* (0.9rho + 0.1 * 50)
        @test b.dissolved ≈ (s.area .* s.H .* poro .- pv_trap) .* 0.8 * 50
        # reachable: the trap a column drains into and those downstream
        r = reachable_capacity(res, b.structural)
        k = findfirst(>(0), res.trap_regions)
        down = downstream_traps(res, res.trap_regions[k])
        @test r[k] ≈ sum(b.structural[findall(in(down), res.traps)])
        @test all(iszero, r[res.trap_regions .== 0])
        # JutulDarcy's CO2 density by default
        c = trapping_capacity(s, res; p, T, res_gas = 0.2, res_water = 0.1,
                              dis_max = 50.0, poro)
        @test all(x -> 100 < x < 1000, c.density)
    end

    @testset "injection case" begin
        mesh, cols = _sine_mesh()
        n = number_of_columns(cols)
        ve = ve_domain(mesh, cols; perm_x = fill(1e-13, n), poro = fill(0.25, n))
        a = (mesh = mesh, columns = cols, ve = ve)
        s = trap_surface(mesh, cols)
        sides, _ = ve_boundary_sides(mesh, cols, ve.domain)
        ns = length(sides)

        # the boundary in order: one loop, consecutive sides a column apart
        loops = ve_boundary_loops(s, sides)
        @test length(loops) == 1 && sort(only(loops)) == 1:ns
        l = only(loops)
        step(i, j) = hypot((sides.edge_midpoint[i] .- sides.edge_midpoint[j])...)
        @test all(step(l[k], l[mod1(k + 1, ns)]) <= 150.0 for k in 1:ns)

        # side types: open and semi-open get conditions, semi-open scaled
        types = [(:closed, :semi_open, :open)[mod1(k, 3)] for k in 1:ns]
        fluid = co2lab_ve_fluid(default_pvt = nothing, residual = (0.1, 0.2))
        yr = 365.2425 * 86400
        c = setup_injection_case(a, [(n ÷ 2, 1e8 / yr)]; fluid, boundary = types,
                                 inj_time = 2yr, inj_steps = 2, mig_time = 10yr,
                                 mig_steps = 1, outside_distance = 1e3)
        @test c.side_types == types && length(c.bc) == count(!=(:closed), types)
        k = findfirst(==(:semi_open), types)
        col = sides.column[k]
        d = hypot(s.cx[col] - sides.edge_midpoint[k][1],
                  s.cy[col] - sides.edge_midpoint[k][2])
        bck = c.bc[findfirst(==(k), findall(!=(:closed), types))]
        @test bck.trans ≈ sides.half_trans[k] * d / (d + 1e3)
        @test c.dt ≈ [yr, yr, 10yr] && c.wells == [n ÷ 2]
        @test_throws ArgumentError setup_injection_case(a, []; fluid, boundary = :leaky)

        # a run handing over each report step; one stopped after the first
        seen = Int[]
        r = simulate_injection_case(c; info_level = -1,
                                    on_report = (k, st) -> (push!(seen, k); true))
        @test seen == [1, 2, 3] && length(r.states) == 3
        @test maximum(r.states[2][:sG]) > 0
        @test isnothing(simulate_injection_case(c; info_level = -1,
                                                on_report = (k, st) -> k < 1))
    end

    @testset "trapping inventory vs co2lab" begin
        c, d = (joinpath(refdir, f) for f in ("inventory_case.mat", "inventory_mrst.mat"))
        if !isfile(c) || !isfile(d) || !isdir(atlas_data_directory())
            @info "Skipping: Atlas data or MRST reference dumps $c, $d missing"
        else
            case = matread(c)
            reps = [_inventory_vs_mrst(case, matread(d)), _inventory_consistency(case)]
            # co2lab's own run of the case, if exported
            run = joinpath(refdir, "inventory_mrst_run.mat")
            isfile(run) && push!(reps, _inventory_run_vs_mrst(case, matread(run)))
            for rep in reps
                passed(rep) || display(rep)
                @test passed(rep)
            end
        end
    end

    @testset "VE fluid vs MRST co2lab" begin
        f = joinpath(refdir, "ve_fluid.mat")
        if !isfile(f)
            @info "Skipping: MRST reference dump $f not found"
        else
            d = matread(f)
            H = vec(d["H"]); g = d["gravity"]
            pvt = (co2_rho_ref = 760.0, wat_rho_ref = 1100.0,
                   co2_rho_pvt = (4.3e-9, 150e5), wat_rho_pvt = (4.3e-10, 100e5),
                   co2_mu_pvt = (2e-9, 150e5))
            # (fluid, caprock rugosity dh): rugosity is a grid property here
            cases = [(co2lab_ve_fluid(; default_pvt = nothing,
                                      residual = (0.0, 0.0), pvt...), 0.0),
                     (co2lab_ve_fluid(; default_pvt = nothing,
                                      residual = (0.27, 0.20), pvt...), 0.0),
                     (co2lab_ve_fluid(; default_pvt = nothing,
                                      residual = (0.27, 0.20), krmax = (0.9, 0.7),
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

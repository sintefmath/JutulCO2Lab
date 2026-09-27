# Comparisons of the trapping inventory with co2lab-ve (roadmap Stage 9 item 2),
# used by the package tests (runtests.jl) and the Stage 9 gate. The states of
# a JutulCO2Lab run (scripts/inventory_case.jl, inventory_case.mat) are given to
# co2lab's massTrappingDistributionVE (scripts/mrst_dump_inventory.m,
# inventory_mrst.mat) and to trapping_inventory here: identical input.

# The model, fluid and traps of the exported case, rebuilt from its settings
# (the inventory needs no wells or forces)
function _inventory_case_model(c)
    case = atlas_ve_domain(c["name"]; coarsening = Int(c["coarsening"]),
                           co2lab_compatible = true, min_thickness = 0.0,
                           rugosity = c["rugosity"])
    fluid = co2lab_ve_fluid(residual = (c["srw"], c["srg"]),
                            co2_rho_ref = c["rho_co2"], wat_rho_ref = c["rho_w"],
                            co2_rho_pvt = Tuple(vec(c["co2_rho_pvt"])),
                            wat_rho_pvt = Tuple(vec(c["wat_rho_pvt"])),
                            dissolution = true)
    model = setup_co2lab_model(case.ve.domain; fluid,
                               rock_compressibility = (reference_pressure = c["pv_ref"],
                                                       compressibility = c["pv_c"]))
    s = trap_surface(case.mesh, case.columns)
    traps = trap_analysis(s; method = :node)
    return (model = model, fluid = fluid, s = s, traps = traps)
end

# Every category at every report step against co2lab's, relative to the CO2
# in place; and the stored inventory of the run itself
function _inventory_vs_mrst(c, d)
    rep = OracleReport("Trapping inventory vs co2lab's massTrappingDistributionVE";
                       mrst = d)
    m = _inventory_case_model(c)
    theirs = d["masses"]
    ours = zeros(size(theirs))
    for k in axes(theirs, 1)
        state = Dict(:Pressure => c["pressure"][:, k], :sG => c["sG"][:, k],
                     :sGmax => c["sGmax"][:, k], :Rs => c["rs"][:, k])
        inv = trapping_inventory(state, m.model, m.fluid, m.s, m.traps)
        ours[k, :] = collect(values(inv))
    end
    scale = maximum(sum(theirs; dims = 2))
    for (j, name) in enumerate(TRAPPING_CATEGORIES)
        check!(rep, "$name", ours[:, j], theirs[:, j];
               tol = Tolerance(0.0, 1e-10 * scale))
    end
    check!(rep, "stored inventory of the run", c["inventory"], ours;
           tol = Tolerance(0.0, 1e-12 * scale))
    return rep
end

# The run's inventory against its own masses: the categories add up to the
# simulator's CO2 in place, and what has left is injected minus in place
function _inventory_consistency(c)
    rep = OracleReport("Trapping inventory: consistency with the simulation")
    total = vec(sum(c["inventory"]; dims = 2))
    scale = maximum(c["injected"])
    check!(rep, "categories sum to the CO2 in place", total, vec(c["co2_in_place"]);
           tol = Tolerance(0.0, 1e-10 * scale))
    left = vec(c["injected"]) .- vec(c["co2_in_place"])
    check!(rep, "exited = injected - in place", vec(c["exited"]), left;
           tol = Tolerance(0.0, 1e-10 * scale))
    return rep
end

# End to end: co2lab's own run of the case on Julia's geometry
# (scripts/mrst_ref_inventory.m, inventory_mrst_run.mat) and its inventory of
# its own states, against the Julia run and its inventory. Tolerances as the
# Stage 4 gate's for runs on identical input (agreed 2026-09-26): 1e-8 —
# relative for pressures, absolute for saturations and rs, relative to the
# injected mass for masses. The injector's BHP is a known well-model
# difference (researcher, 2026-09-28): JutulDarcy converts the injected
# CO2's surface rate to reservoir volume at the wellbore pressure, co2lab at
# the cell pressure, so with compressible CO2 the drawdown differs by c_G
# times the drawdown; the check is that the difference stays within that.
function _inventory_run_vs_mrst(c, r; tol = 1e-8)
    rep = OracleReport("Inventory case: Julia run vs co2lab's run, end to end"; mrst = r)
    rel, abs_ = Tolerance(tol, 0.0), Tolerance(0.0, tol)
    check!(rep, "pressure", c["pressure"], r["pressure"]; tol = rel)
    for k in ("sG", "sGmax", "rs")
        check!(rep, k, c[k], r[k]; tol = abs_)
    end
    on = findall(>(0), vec(c["injected"]) .- [0.0; vec(c["injected"])[1:end-1]])
    _bhp_vs_mrst!(rep, c, r, on)

    scale = Tolerance(0.0, tol * maximum(c["injected"]))
    for (j, name) in enumerate(TRAPPING_CATEGORIES)
        check!(rep, "inventory: $name", c["inventory"][:, j], r["masses"][:, j];
               tol = scale)
    end
    check!(rep, "injected", vec(c["injected"]), vec(r["injected"]); tol = scale)
    check!(rep, "exited", vec(c["exited"]), vec(r["exited"]); tol = scale)
    return rep
end

# The injector's BHP against co2lab's: the difference, relative to c_G times
# the drawdown (BHP minus the well column's pressure), must not exceed 1.
# Modifies: `rep` (a check added).
function _bhp_vs_mrst!(rep, c, r, on)
    col = Int(c["well_column"])
    bhp, bhp_mrst = vec(c["bhp"])[on], vec(r["bhp"])[on]
    drawdown = bhp .- c["pressure"][col, on]
    diff = abs.(bhp .- bhp_mrst)
    ratio = diff ./ (c["co2_rho_pvt"][1] .* drawdown .^ 2)
    rel = round(maximum(diff ./ bhp_mrst); sigdigits = 2)
    check!(rep, "injector BHP difference / (c_G × drawdown²) (well model; max " *
           "rel $rel)", maximum(ratio), 0.0; tol = Tolerance(0.0, 1.0))
    return rep
end

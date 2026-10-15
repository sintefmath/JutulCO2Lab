# The capacity estimate against co2lab-demo's exploreCapacity, used by the
# package tests (runtests.jl). Input: the MRST export capacity.mat.

using SparseArrays

# A case of the dump as trapping_capacity reads it: the columns (depth,
# height, area) and the trap analysis
function _capacity_case(e)
    v(x) = vec(Float64.(x))
    iv(x) = vec(round.(Int, x))
    s = (zc = v(e["z"]), H = v(e["H"]), area = v(e["area"]))
    adj = sparse(e["trap_adj"] isa Number ? fill(e["trap_adj"], 1, 1) : e["trap_adj"])
    res = (traps = iv(e["traps"]), trap_z = v(e["trap_z"]),
           trap_regions = iv(e["trap_regions"]), trap_adj = adj .!= 0)
    return (s, res)
end

# Each case's capacities and reachable capacity, on co2lab's columns, traps,
# pressure, temperature and CO2 density (identical input)
function _capacity_vs_mrst(d)
    rep = OracleReport("capacity estimate vs co2lab-demo's exploreCapacity " *
                       "(identical input)")
    v(x) = vec(Float64.(x))
    tol = Tolerance(1e-12, 1e-6)   # kg
    for (k, e) in enumerate(vec(d["cases"]))
        s, res = _capacity_case(e)
        label = e["name"] * (isodd(k) ? " (defaults)" : " (other values)")
        # co2lab's dis_max is in reference volumes, ours in kg per m3 brine
        cap = trapping_capacity(s, res; p = v(e["p"]), T = v(e["t"]),
                                res_gas = e["res_sat_co2"], res_water = e["res_sat_wat"],
                                dis_max = e["dis_max"] * d["rhoCref"],
                                poro = v(e["poro"]), co2_density = v(e["rho"]))
        check!(rep, "$label: structural", cap.structural, v(e["strap"]); tol)
        check!(rep, "$label: residual", cap.residual, v(e["btrap_res"]); tol)
        check!(rep, "$label: dissolved", cap.dissolved, v(e["btrap_dis"]); tol)
        check!(rep, "$label: reachable structural",
               reachable_capacity(res, cap.structural), v(e["cum_cap"]); tol)
    end
    return rep
end

# Comparisons of the capillary-fringe vertical model with co2lab's
# 'P-scaled table' fluid (roadmap Stage 6 item 1 and its pre-gate), used by
# the package tests (runtests.jl) and the Stage 6 gate. Reference:
# scripts/mrst_dump_fringe_fluid.m (fringe_fluid.mat): co2lab's sampled table
# and its krG, krW, pcWG and heights on sampled states, per case.

# The fluid of a dumped case: the PVT of the dump script, and the case's
# fine-scale curves (the slope cases' (C, alpha), Johansen's functions, the
# Leverett cases' J-curve)
function _fringe_case_fluid(c)
    name, rw, rg = c["name"], c["res_water"], c["res_gas"]
    if startswith(name, "slope")
        curves = (invPc3D = (0.4, 0.5), kr3D = s -> s^2)
    elseif startswith(name, "Leverett")
        curves = (invPc3D = J -> (1 - rw) * (0.5 / (J + 0.5))^2 + rw, kr3D = s -> s^2,
                  leverett = true)
    else
        pe = 5e3
        curves = (invPc3D = pc -> (1 - rw) * (pe / max(pc, pe))^2 + rw,
                  kr3D = s -> max((s - rg) / (1 - rg), 0.0)^2)
    end
    return co2lab_ve_fluid(; vertical = CapillaryFringeTable(; curves...),
                           residual = (rw, rg), co2_rho_ref = 760.0,
                           wat_rho_ref = 1100.0, co2_rho_pvt = (4.3e-9, 150e5),
                           wat_rho_pvt = (4.3e-10, 100e5), co2_mu_ref = 6e-5,
                           wat_mu_ref = 8e-4)
end

# Our table and our krG, krW, pcWG, heights for every sample of case `c`
function _fringe_case_values(c, H)
    fluid = _fringe_case_fluid(c)
    domain = Dict(:column_height => H, :porosity => vec(c["poro"]),
                  :permeability => vec(c["perm"]))
    m = JutulCO2Lab._vertical_model(fluid.vertical, domain, fluid)
    rw, rg, g = fluid.res_water, fluid.res_gas, Jutul.gravity_constant
    out = Dict(k => zeros(size(c["sg"])) for k in ("krG", "krW", "pcWG", "h", "hmax"))
    for s in axes(c["sg"], 2), i in eachindex(H)
        sg, smax, p = c["sg"][i, s], c["sgmax"][i, s], c["p"][i, s]
        drho = JutulCO2Lab._density_difference(fluid, p, 0.0)
        h, hmax = JutulCO2Lab._column_heights(m, i, sg, smax, H[i], drho, rw, rg)
        krW, krG = JutulCO2Lab._column_relperm(m, i, sg, smax, h, hmax, H[i], drho, 0.0,
                                               fluid.krw, fluid.krg, rw, rg)
        out["krG"][i, s], out["krW"][i, s] = krG, krW
        out["pcWG"][i, s], out["h"][i, s], out["hmax"][i, s] = drho * g * h, h, hmax
    end
    return (m = m, values = out)
end

# Table and values of every case against co2lab's ('P-scaled table') or
# co2lab's corrected by scripts/mrst_pk_patch ('P-K-scaled table'). Plumes truncated at the
# column bottom are compared at `tol_trunc`: co2lab stops its Newton solve for
# them at a residual of 1e-8 times the largest Δρ g H of all columns (about
# 1e-7 relative in Pc), where ours iterates to round-off.
function _fringe_fluid_vs_mrst(d; tol = 1e-10, tol_trunc = 1e-6)
    rep = OracleReport("Fringe fluid ('P-scaled table') vs co2lab's makeVEFluid";
                       mrst = d)
    H = vec(d["H"])
    for k in 1:Int(d["ncases"])
        c = d["case$k"]
        ours = _fringe_case_values(c, H)
        label = c["name"]
        # Leverett: each column's factor σ√(φ/k)
        haskey(c, "kscale") && check!(rep, "$label: kscale", ours.m.kscale,
                                      vec(c["kscale"]); tol = Tolerance(1e-14, 0.0))
        for (field, ref) in (("p", "table_p"), ("SP", "table_SP"), ("krP", "table_krP"))
            check!(rep, "$label: table $field", getfield(ours.m, Symbol(field)),
                   vec(c[ref]); tol = Tolerance(1e-12, 1e-12 * maximum(c[ref])))
        end
        # co2lab's truncated plumes: heights below the column (h_max: its own)
        for q in ("pcWG", "krG", "krW", "h", "hmax")
            trunc = (q == "hmax" ? c["hmax"] : c["h"]) .> H
            for (sel, t, tag) in ((.!trunc, tol, ""),
                                  (trunc, tol_trunc, " (truncated plumes)"))
                ref = c[q][sel]
                isempty(ref) && continue
                check!(rep, "$label: $q$tag", ours.values[q][sel], ref;
                       tol = Tolerance(t, t * max(maximum(abs, ref), 1e-12)))
            end
        end
    end
    return rep
end

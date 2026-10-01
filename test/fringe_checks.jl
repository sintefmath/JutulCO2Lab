# Comparisons of the capillary-fringe vertical model with co2lab's
# 'P-scaled table' and 'P-K-scaled table' fluids (roadmap Stage 6 item 1 and
# its pre-gate), with JutulCO2Lab's patch of co2lab's fringe tables
# (scripts/mrst_fringe_patch: adaptive Hermite tables, plume-tip
# regularisation, the corrected 'P-K-scaled table'), used by the package tests
# (runtests.jl) and the Stage 6 gate. Reference: scripts/mrst_dump_fringe_fluid.m
# (fringe_fluid.mat): co2lab's table and its krG, krW, pcWG and heights on
# sampled states, per case.

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

# Table and values of every case against the patched co2lab's. Both sides
# solve for Pc to round-off, plumes reaching the bottom included.
function _fringe_fluid_vs_mrst(d; tol = 1e-10, tol_trunc = tol)
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
        check!(rep, "$label: table nodes", length(ours.m.p), length(vec(c["table_p"]));
               tol = Tolerance(0.0, 0.0))
        length(ours.m.p) == length(vec(c["table_p"])) || continue
        for (field, ref) in (("p", "table_p"), ("SP", "table_SP"), ("krP", "table_krP"),
                             ("s", "table_dSP"), ("kr", "table_dkrP"))
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

# capfringeFigures.m, capillary-fringe half (mrst_dump_capfringe_cp.m,
# capfringe_cp.mat): the 'P-scaled table' curves over sG with and without
# hysteresis, the heights and the reconstructed profile for S = 0.2,
# Smax = 0.3, on the example's 50 m column. co2lab's table uses the reference
# densities; the curves the densities at the example's pressure, reproduced by
# linear PVT through them. Plumes reaching the bottom at `tol_trunc`
# (co2lab's Newton accuracy there). A full column (S = 1 − srw) has no
# determined Pc — any large enough one fills it — so there the check is that
# both Pc give back S.
function _capfringe_cp_vs_mrst(d; tol = 1e-10, tol_trunc = tol)
    J = JutulCO2Lab
    mesh = ve_mesh_from_mrst(d["G3"])
    cols, _ = ve_columns(mesh)
    n = Int(d["G3"]["cells"]["num"])
    ve = ve_domain(mesh, cols; perm_x = fill(1e-13, n), poro = fill(0.3, n))
    H, rw, rg, p = ve.geometry.H[1], d["srw"], d["src"], d["pressure"]
    # linear PVT through co2lab's densities at the example's pressure
    lin(rho, ref) = ((rho / ref - 1) / p, 0.0)
    fringe = CapillaryFringeTable(; invPc3D = (d["C"], d["alpha"]),
                                  kr3D = s -> max((s - rg) / (1 - rg), 0.0)^2)
    fl = co2lab_ve_fluid(; vertical = fringe, residual = (rw, rg),
                         wat_rho_ref = d["rhoWS"], co2_rho_ref = d["rhoGS"],
                         wat_rho_pvt = lin(d["rhoW"], d["rhoWS"]),
                         co2_rho_pvt = lin(d["rhoG"], d["rhoGS"]))
    m = J._vertical_model(fl.vertical, ve.domain, fl)
    drho = J._density_difference(fl, p, 0.0)
    g = Jutul.gravity_constant

    function curves(s, smax)
        hh = J._column_heights.(Ref(m), 1, s, smax, H, drho, rw, rg)
        kr = [J._column_relperm(m, 1, s[i], smax[i], hh[i]..., H, drho, 0.0, fl.krw,
                                fl.krg, rw, rg) for i in eachindex(s)]
        return first.(kr), last.(kr), drho * g .* first.(hh), first.(hh)
    end

    rep = OracleReport("capfringeFigures.m, capillary fringe, vs MRST"; mrst = d)
    s = vec(d["s"])
    full = s .>= 1 - rw - 1e-12
    hyst = (" with hysteresis", vec(d["smax_hyst"]), "_hyst")
    for (label, smax, sfx) in (("", s, ""), hyst)
        krw, krg, pc, h = curves(s, smax)
        trunc = (h .> H) .& .!full
        rest = .!(h .> H) .& .!full
        groups = ((rest, tol, ""), (trunc, tol_trunc, ", plumes reaching the bottom"))
        for (sel, t, tag) in groups
            any(sel) || continue
            for (q, ours) in (("krW", krw), ("krG", krg), ("pcWG", pc))
                ref = vec(d[q * sfx])[sel]
                check!(rep, "$q$label$tag", ours[sel], ref;
                       tol = Tolerance(t, t * maximum(abs, ref; init = 1e-12)))
            end
        end
        # the full column: the saturation each Pc gives back
        dP = drho * g * H
        sback(q) = (J._SP(m, q) - J._SP(m, q - dP)) / dP
        check!(rep, "full column$label: S from Pc, ours and co2lab's",
               [sback.(pc[full]); sback.(vec(d["pcWG" * sfx])[full])],
               [s[full]; s[full]]; tol = Tolerance(0.0, 1e-8))
    end

    h, hmax = J._column_heights(m, 1, d["S"], d["Smax"], H, drho, rw, rg)
    check!(rep, "profile: heights", [h, hmax], [d["h"], d["h_max"]];
           tol = Tolerance(tol, 0.0))
    st = Dict(:Pressure => [p], :sG => [d["S"]], :sGmax => [d["Smax"]])
    o = ve_reconstruct(st, ve, cols, fl)
    check!(rep, "profile: sG", o.sG, vec(d["s3D"]); tol = Tolerance(0.0, tol))
    check!(rep, "profile: sGmax", o.sGmax, vec(d["smax3D"]); tol = Tolerance(0.0, tol))
    return rep
end

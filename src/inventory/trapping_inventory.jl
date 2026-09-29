# The trapping inventory of a VE simulation (roadmap Stage 9 item 2):
# co2lab-ve's massTrappingDistributionVE and postprocessStates — where the CO2
# is, by trapping mechanism, per state and over a run.
#
# Per column, the plume (height h) holds mobile CO2 above the residual zone
# (down to h_max). Both are split at the column's trap height zt (from the
# caprock down to its trap's spill depth, from trap_analysis) into structural
# and free parts, and at the subscale height dh (the caprock's rugosity) into
# subscale-trapped parts; the dissolved CO2 is a category of its own.

"Names of the inventory categories, in co2lab's reporting order."
const TRAPPING_CATEGORIES = (:dissolved, :structural_residual, :residual,
                             :residual_in_plume, :structural_subscale,
                             :structural_plume, :free_plume)

# ------------------------------ Column properties -----------------------------

# What the inventory needs per column besides the state: static pore volume,
# height, rugosity (the grid's, if it has one), trap height (from the trap
# analysis `res` on the trap surface `s`), and the pressure dependence of
# pore volume and densities
function _inventory_columns(model, fluid::CO2LabVEFluid, s, res)
    rmodel = JutulDarcy.reservoir_model(model)
    domain = JutulDarcy.reservoir_domain(rmodel)
    H = domain[:column_height]
    dh = haskey(domain, :rugosity, Cells()) ? domain[:rugosity] : zeros(length(H))
    length(s.zc) == length(H) ||
        throw(ArgumentError("the trap surface has $(length(s.zc)) columns, the " *
                            "model $(length(H))"))

    # the height below the caprock that lies in a trap (co2lab's zt)
    spill(r) = r == 0 ? 0.0 : res.trap_z[r]
    zt = [min(max(spill(r) - z, 0.0), h) for (r, z, h) in zip(res.trap_regions, s.zc, H)]
    return (pv = JutulDarcy.pore_volume(domain), H = H, dh = dh, zt = zt,
            pv_mult = _pore_volume_multiplier(rmodel))
end

# The pore volume's pressure dependence as a function of pressure: co2lab's
# linear form if set (set_co2lab_rock_compressibility!), else none
function _pore_volume_multiplier(rmodel)
    v = get(rmodel.secondary_variables, :FluidVolume, nothing)
    isnothing(v) && return p -> 1.0
    v isa LinearPoreVolume && return p -> 1 + v.c * (p - v.p_ref)
    throw(ArgumentError("pore volume model $(typeof(v)) not supported by the " *
                        "inventory; use co2lab's rock compressibility"))
end

# ---------------------------------- One state ---------------------------------

"""
    trapping_inventory(state, model, fluid, s, res) -> NamedTuple

Where the CO2 of a VE simulation state is, by trapping mechanism — co2lab's
`massTrappingDistributionVE` for a sharp interface. `state` is a reservoir
state of a run of `model` ([`setup_co2lab_model`](@ref) with `fluid`); `s`
and `res` are the trap surface ([`trap_surface`](@ref)) and trap analysis
([`trap_analysis`](@ref)) of the model's columns. The subscale height is the
grid's rugosity (`:rugosity` of the VE domain), zero if it has none.

Returns a new NamedTuple of masses (kg), in co2lab's reporting order
([`TRAPPING_CATEGORIES`](@ref)):
- `dissolved`: dissolved in the brine;
- `structural_residual`: residual CO2 within traps;
- `residual`: residual CO2 outside traps, below the plume;
- `residual_in_plume`: the part of the plume outside traps that will stay
  behind as residual CO2 when the plume moves on;
- `structural_subscale`: in the caprock's subscale relief (rugosity);
- `structural_plume`: mobile CO2 within traps;
- `free_plume`: mobile CO2 outside traps.
Each is at least zero, as in co2lab. The arguments are only read.
"""
function trapping_inventory(state, model, fluid::CO2LabVEFluid, s, res)
    return _trapping_inventory(state, fluid, _inventory_columns(model, fluid, s, res))
end

# The inventory of a state with the column properties `col` already assembled
function _trapping_inventory(state, fluid, col)
    rw, rg = fluid.res_water, fluid.res_gas
    (cw, pw), (cg, pg) = fluid.rho_pvt
    p, S, Smax = state[:Pressure], state[:sG], state[:sGmax]
    # no Rs in the state: a model without dissolution
    rs = haskey(state, :Rs) ? state[:Rs] : zeros(length(p))

    # pore volume and CO2 density at the state's pressure
    pv = col.pv .* col.pv_mult.(p)
    rho_ref = fluid.rho_ref[2]
    mfac = pv .* rho_ref .* (1 .+ cg .* (p .- pg))
    # brine expansion factor: rs is surface CO2 per surface brine volume
    bW = 1 .+ cw .* (p .- pw)

    hh = [_plume_heights(S[i], Smax[i], col.H[i], rw, rg) for i in eachindex(S)]
    # co2lab adds eps to avoid division by zero
    h, h_max = first.(hh) .+ eps(), last.(hh) .+ eps()
    m = _category_saturations(S, Smax, h, h_max, col.zt, col.dh, rw, rg)

    # residual below the plume; the rest of the free immobile part is in the plume
    residual = sum(mfac .* m.res)
    masses = (dissolved = sum(rho_ref .* bW .* pv .* rs .* (1 .- S)),
              structural_residual = sum(mfac .* m.immob_struct),
              residual = residual,
              residual_in_plume = sum(mfac .* m.immob_free) - residual,
              structural_subscale = sum(mfac .* (m.immob_sub .+ m.mobile_sub)),
              structural_plume = sum(mfac .* m.mobile_struct),
              free_plume = sum(mfac .* m.mobile_free))
    return map(x -> max(x, 0.0), masses)
end

# co2lab's split of each column's CO2 saturation: immobile (residual) and
# mobile parts, each in subscale traps (dh), in structural traps (zt) and free;
# `res` the residual CO2 outside traps below the plume
function _category_saturations(S, Smax, h, h_max, zt, dh, rw, rg)
    # the saturation left when the mobile part is gone
    S_immob = Smax .* rg ./ (1 - rw)
    S_mob = S .- S_immob
    immob_sub = S_immob .* min.(dh, h_max) ./ h_max
    # each part scaled by the fraction of its height in the zone
    mobile_sub = S_mob .* min.(dh, h) ./ h
    immob_struct = S_immob .* max.(min.(h_max, zt) .- dh, 0) ./ h_max
    mobile_struct = S_mob .* max.(min.(h, zt) .- dh, 0) ./ h
    return (immob_sub = immob_sub, mobile_sub = mobile_sub,
            immob_struct = immob_struct, mobile_struct = mobile_struct,
            immob_free = S_immob .- immob_struct .- immob_sub,
            mobile_free = S_mob .- mobile_struct .- mobile_sub,
            res = S_immob .* max.(h_max .- max.(h, zt), 0) ./ h_max)
end

# ---------------------------------- A run -------------------------------------

"""
    trapping_history(result, model, fluid, s, res; initial = 0.0) -> NamedTuple

The trapping inventory ([`trapping_inventory`](@ref)) of every report step of
a run (co2lab's `postprocessStates`), with the CO2 that has left the domain:
`result` is `simulate_reservoir`'s result (a `ReservoirSimResult`), the
other arguments as for [`trapping_inventory`](@ref); `initial` the CO2 mass
in place at the start.

Returns the new NamedTuple `(t, inventory, injected, exited)`: the report
times (s), one inventory per step, the CO2 mass injected up to each step
(net, from the wells' CO2 mass rates) and the mass that has left the domain
(`initial + injected` minus the inventory's total). The arguments are only
read.
"""
function trapping_history(result, model, fluid::CO2LabVEFluid, s, res;
                          initial::Real = 0.0)
    col = _inventory_columns(model, fluid, s, res)
    inventory = [_trapping_inventory(st, fluid, col) for st in result.states]
    t = result.time
    dt = diff([0.0; t])
    # the wells' CO2 (vapour) mass rates, positive into the domain
    rate = zeros(length(t))
    for w in values(result.wells.wells)
        rate .+= w[:Vapor_mass_rate]
    end
    injected = cumsum(rate .* dt)
    # what is neither in place nor unaccounted for has left through the boundary
    exited = initial .+ injected .- sum.(values.(inventory))
    return (t = t, inventory = inventory, injected = injected, exited = exited)
end

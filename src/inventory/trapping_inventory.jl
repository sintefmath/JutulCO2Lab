# The trapping inventory of a VE simulation: where the CO2 is, by trapping
# mechanism, per state and over a run (co2lab-ve's massTrappingDistributionVE
# and postprocessStates).
#
# Per column, the plume (height h) holds mobile CO2 above the residual zone
# (down to h_max). Both are split at the column's trap height zt (caprock down
# to its trap's spill depth, from trap_analysis) into structural and free
# parts, and at the subscale height dh (the caprock's rugosity) into
# subscale-trapped parts. Dissolved CO2 is a category of its own.

"Names of the inventory categories, in co2lab's reporting order."
const TRAPPING_CATEGORIES = (:dissolved, :structural_residual, :residual,
                             :residual_in_plume, :structural_subscale,
                             :structural_plume, :free_plume)

# ------------------------------ Column properties -----------------------------

# Per-column inputs of the inventory besides the state: static pore volume,
# height, rugosity (the grid's, if any), trap height (from trap analyses `res`
# on trap surfaces `s`: one of each, or vectors, one per sheet), pore-volume
# pressure dependence, temperature (tabulated properties) and the vertical
# model
function _inventory_columns(model, fluid::CO2LabVEFluid, s, res)
    rmodel = JutulDarcy.reservoir_model(model)
    domain = JutulDarcy.reservoir_domain(rmodel)
    H = domain[:column_height]
    dh = haskey(domain, :rugosity, Cells()) ? domain[:rugosity] : zeros(length(H))
    zt = s isa AbstractVector ? _trap_heights(H, s, res) : _trap_heights(H, [s], [res])
    T = haskey(domain, :temperature, Cells()) ? domain[:temperature] :
        fill(Jutul.default_value(nothing, JutulDarcy.Temperature()), length(H))
    return (pv = JutulDarcy.pore_volume(domain), H = H, dh = dh, zt = zt, T = T,
            pv_mult = _pore_volume_multiplier(rmodel),
            vertical = _vertical_model(fluid.vertical, domain, fluid))
end

# Per column, the height below its top that lies in a trap (co2lab's zt),
# from each surface's polygons (its `blocks`; without them, polygon i is
# column i and the surface must cover every column) and trap regions;
# 0 for columns on no surface (fine cells, sheets not given)
function _trap_heights(H, surfaces, results)
    length(surfaces) == length(results) ||
        throw(ArgumentError("one trap analysis per trap surface needed"))
    zt = zeros(length(H))
    for (s, res) in zip(surfaces, results)
        blocks = hasproperty(s, :blocks) ? s.blocks : eachindex(s.zc)
        (hasproperty(s, :blocks) ? maximum(blocks; init = 0) <= length(H) :
                                   length(s.zc) == length(H)) ||
            throw(ArgumentError("the trap surface does not fit the model's " *
                                "$(length(H)) columns"))
        spill(r) = r == 0 ? 0.0 : res.trap_z[r]
        for (i, c) in enumerate(blocks)
            zt[c] = min(max(spill(res.trap_regions[i]) - s.zc[i], 0.0), H[c])
        end
    end
    return zt
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
    trapping_inventory(state, model, fluid, s, res; co2lab_compatible = false)
        -> NamedTuple

Where the CO2 of one VE simulation state is, by trapping mechanism (co2lab's
`massTrappingDistributionVE`). [`trapping_history`](@ref) does this for a
whole run.

# Arguments

- `state`: a reservoir state of a run of `model`.
- `model`, `fluid`: the model ([`setup_co2lab_model`](@ref)) and its fluid.
- `s`, `res`: a trap surface ([`trap_surface`](@ref)) and its
  [`trap_analysis`](@ref); or vectors of them, one per sheet
  ([`ve_sheets`](@ref)) of a grid with several blocks per pillar.
- `co2lab_compatible`: for a capillary fringe, split by height as co2lab
  does (see Notes).

# Returns

A new NamedTuple of CO2 masses (kg), in co2lab's reporting order
([`TRAPPING_CATEGORIES`](@ref)):
- `dissolved`: dissolved in the brine;
- `structural_residual`: residual CO2 within traps;
- `residual`: residual CO2 outside traps, below the plume;
- `residual_in_plume`: the part of the plume outside traps that will stay
  behind as residual CO2 when the plume moves on;
- `structural_subscale`: in the caprock's subscale relief (rugosity);
- `structural_plume`: mobile CO2 within traps;
- `free_plume`: mobile CO2 outside traps.

Each is at least zero, as in co2lab.

# Notes

- Columns on none of the surfaces (fine cells, sheets not given) have no
  structural trapping.
- The subscale height is the VE domain's `:rugosity`, zero if it has none.
- Each column's immobile (residual) and mobile CO2 is split between the
  zones under the caprock (subscale relief, structural trap, below) by
  where it lies in the column:
  - sharp interface: the saturation is even over the plume, so by height,
    as co2lab;
  - capillary fringe: by the fringe's saturation profile (mobile part by
    the current one, immobile part by the maximum one). co2lab splits by
    height here too, which is what `co2lab_compatible = true` does.
"""
function trapping_inventory(state, model, fluid::CO2LabVEFluid, s, res;
                            co2lab_compatible::Bool = false)
    col = _inventory_columns(model, fluid, s, res)
    return _trapping_inventory(state, fluid, col; co2lab_compatible)
end

# The inventory of a state with the column properties `col` already assembled
function _trapping_inventory(state, fluid, col; co2lab_compatible::Bool = false)
    rw, rg = fluid.res_water, fluid.res_gas
    p, S, Smax = state[:Pressure], state[:sG], state[:sGmax]
    # no Rs in the state: a model without dissolution
    rs = haskey(state, :Rs) ? state[:Rs] : zeros(length(p))

    # pore volume and CO2 density at the state's pressure (and temperature)
    pv = col.pv .* col.pv_mult.(p)
    rho_ref = fluid.rho_ref[2]
    b = [_shrinkage(fluid, p[i], col.T[i]) for i in eachindex(p)]
    mfac = pv .* rho_ref .* last.(b)
    # brine expansion factor: rs is surface CO2 per surface brine volume
    bW = first.(b)

    # heights from the vertical model, at each column's density difference
    m = col.vertical
    drho = [_density_difference(fluid, p[i], rs[i], col.T[i]) for i in eachindex(p)]
    hh = [_column_heights(m, i, S[i], Smax[i], col.H[i], drho[i], rw, rg)
          for i in eachindex(S)]
    c = _category_saturations(m, S, Smax, first.(hh), last.(hh), col, drho, rw, rg;
                              co2lab_compatible)

    # residual below the plume; the rest of the free immobile part is in the plume
    residual = sum(mfac .* c.res)
    masses = (dissolved = sum(rho_ref .* bW .* pv .* rs .* (1 .- S)),
              structural_residual = sum(mfac .* c.immob_struct),
              residual = residual,
              residual_in_plume = sum(mfac .* c.immob_free) - residual,
              structural_subscale = sum(mfac .* (c.immob_sub .+ c.mobile_sub)),
              structural_plume = sum(mfac .* c.mobile_struct),
              free_plume = sum(mfac .* c.mobile_free))
    return map(x -> max(x, 0.0), masses)
end

"""
    _category_saturations(m, S, Smax, h, h_max, col, drho, rw, rg;
                          co2lab_compatible = false) -> NamedTuple

Each column's CO2 saturation split by the vertical model `m` into immobile
(residual) and mobile parts, each in subscale traps (depth `col.dh`), in
structural traps (`col.zt`) and free. For the sharp interface this is
co2lab's `massTrappingDistributionVE`.

`h`, `h_max` are the model's heights (uncapped), `drho` the density
differences, `col` the inventory's column properties.

Returns a new NamedTuple of vectors, one entry per column: `immob_sub`,
`mobile_sub`, `immob_struct`, `mobile_struct`, `immob_free`, `mobile_free`,
and `res`, the residual CO2 outside traps below the plume.
"""
function _category_saturations(::SharpInterface, S, Smax, h, h_max, col, drho, rw, rg;
                               co2lab_compatible::Bool = false)
    return _height_split(S, Smax, h, h_max, col, rw, rg)
end

# co2lab's split by height: the immobile part spread evenly down to h_max, the
# mobile part down to h. Returns a new NamedTuple of vectors.
function _height_split(S, Smax, h, h_max, col, rw, rg)
    zt, dh = col.zt, col.dh
    # capped at the column height (a fringe's can exceed it), then co2lab adds
    # eps to avoid division by zero
    h = min.(h, col.H) .+ eps()
    h_max = min.(h_max, col.H) .+ eps()
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
    trapping_history(result, model, fluid, s, res; initial = 0.0,
                     co2lab_compatible = false) -> NamedTuple

The trapping inventory ([`trapping_inventory`](@ref)) at every report step
of a run, with the CO2 that has left the domain (co2lab's
`postprocessStates`).

# Arguments

- `result`: `simulate_reservoir`'s result (a `ReservoirSimResult`).
- `model`, `fluid`, `s`, `res`, `co2lab_compatible`: as for
  [`trapping_inventory`](@ref).
- `initial`: the CO2 mass (kg) in place at the start.

# Returns

A new NamedTuple:
- `t`: the report times (s);
- `inventory`: one [`trapping_inventory`](@ref) per step;
- `injected`: net CO2 mass (kg) injected up to each step, from the wells'
  CO2 mass rates;
- `exited`: CO2 mass (kg) that has left the domain, `initial + injected`
  minus the inventory's total.
"""
function trapping_history(result, model, fluid::CO2LabVEFluid, s, res;
                          initial::Real = 0.0, co2lab_compatible::Bool = false)
    col = _inventory_columns(model, fluid, s, res)
    inventory = [_trapping_inventory(st, fluid, col; co2lab_compatible)
                 for st in result.states]
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

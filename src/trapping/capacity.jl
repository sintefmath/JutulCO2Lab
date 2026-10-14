# A formation's CO2 storage capacity, column by column: co2lab-demo's
# exploreCapacity estimate (its recompute_trapcap and compute_cumul_trapcap),
# an upper bound with every trap filled and all the brine saturated.

# JutulDarcy's CO2 density table, made at first use (building it takes a while)
const _CO2_DENSITY = Ref{Any}(nothing)

# The default CO2 density at (p, T): JutulDarcy's tables
function _default_co2_density()
    if isnothing(_CO2_DENSITY[])
        _CO2_DENSITY[] = jutuldarcy_property_tables().density[2]
    end
    return _CO2_DENSITY[]
end

"""
    trapping_capacity(s, res; p, T, res_gas, res_water, dis_max,
                      poro = ones(ncol), co2_density = nothing) -> NamedTuple

The CO2 a formation can store, per column, as co2lab-demo's
`exploreCapacity` estimates it: an upper bound, with every structural trap
filled to its spill depth and all the brine saturated with dissolved CO2.

# Arguments

- `s`: the [`trap_surface`](@ref).
- `res`: its [`trap_analysis`](@ref).
- `p`, `T`: caprock pressure (Pa) and temperature (K) per column, for the
  CO2 density.
- `res_gas`, `res_water`: residual CO2 and brine saturations.
- `dis_max`: dissolved CO2 mass per brine volume at saturation (kg/m³;
  co2lab gives it in reference volumes: multiply by its reference density
  of 760 kg/m³).
- `poro`: porosity per column, or one value for all; default 1.
- `co2_density`: a function `(p, T) -> kg/m³` (e.g. a [`PhaseTable`](@ref)),
  a vector per column, or `nothing` for JutulDarcy's tables (the CO2
  density of [`jutuldarcy_property_tables`](@ref)).

# Returns

A new NamedTuple, one value per column:
- `structural`: CO2 mass (kg) in the trap part: free CO2 in the pore volume
  less the residual brine (`1 − res_water`), plus CO2 dissolved in that
  residual brine;
- `residual`: residual CO2 mass (kg) below it (`res_gas` of the pore volume);
- `dissolved`: CO2 mass (kg) dissolved below it, in the brine
  (`1 − res_gas` of the pore volume);
- `trap_height`: height of the trap part (m);
- `density`: the CO2 density used (kg/m³).

Each column is split at its trap's spill depth: the trap part runs from the
caprock down to the spill depth (at most the column height `H`), the rest
lies below.
"""
function trapping_capacity(s, res; p, T, res_gas::Real, res_water::Real,
                           dis_max::Real, poro = ones(length(s.zc)),
                           co2_density = nothing)
    n = length(s.zc)
    rho = co2_density isa AbstractVector ? collect(Float64, co2_density) :
          map(something(co2_density, _default_co2_density()), p, T)
    poro = poro isa Real ? fill(Float64(poro), n) : poro
    # each column's part in a trap: down to the trap's spill depth, within H
    h = [t == 0 ? 0.0 : min(res.trap_z[t] - s.zc[k], s.H[k])
         for (k, t) in enumerate(res.traps)]
    pv_trap = s.area .* h .* poro
    pv_below = s.area .* (s.H .- h) .* poro
    return (structural = pv_trap .* ((1 - res_water) .* rho .+ res_water * dis_max),
            residual = pv_below .* res_gas .* rho,
            dissolved = pv_below .* (1 - res_gas) * dis_max,
            trap_height = h, density = rho)
end

"""
    reachable_capacity(res, structural) -> Vector{Float64}

Per column, the structural capacity that CO2 injected there could fill
(co2lab-demo's reachable structural capacity): that of the trap the column
drains into and of every trap downstream of it.

# Arguments

- `res`: the [`trap_analysis`](@ref) result.
- `structural`: structural capacity per column (kg), e.g. from
  [`trapping_capacity`](@ref).

# Returns

A new vector, one value per column (kg); 0 for columns draining out of the
domain.
"""
function reachable_capacity(res, structural)
    ntrap = length(res.trap_z)
    # each trap's structural capacity, then summed down each spill path
    per_trap = zeros(ntrap)
    for (k, t) in enumerate(res.traps)
        t > 0 && (per_trap[t] += structural[k])
    end
    reach = [sum(per_trap[downstream_traps(res, t)]) for t in 1:ntrap]
    return [t == 0 ? 0.0 : reach[t] for t in res.trap_regions]
end

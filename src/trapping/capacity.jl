# A formation's CO2 storage capacity, column by column (roadmap Stage 11
# item 2): co2lab-demo's exploreCapacity estimate (its recompute_trapcap and
# compute_cumul_trapcap) — an upper bound, every trap filled and all the brine
# saturated.

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
`s` is the [`trap_surface`](@ref), `res` its [`trap_analysis`](@ref).

Each column is split at the traps' spill depth: its part in a trap (from the
caprock down to the trap's spill depth, at most its height `H`) and the part
below.
- `structural`: in the trap part, free CO2 in the pore volume less the
  residual brine (`1 − res_water`), plus CO2 dissolved in that residual brine;
- `residual`: below, residual CO2 (`res_gas` of the pore volume);
- `dissolved`: below, CO2 dissolved in the rest of the pore volume, the brine
  (`1 − res_gas`).

`p` and `T`: the caprock pressure (Pa) and temperature (K) per column, at
which the CO2 density is taken. `dis_max`: dissolved CO2 mass per brine
volume when saturated (kg/m³; co2lab gives it in reference volumes, times
its reference density of 760 kg/m³). `poro`: porosity per column.
`co2_density`: a function `(p, T) -> kg/m³` (e.g. a [`PhaseTable`](@ref)), a
vector per column, or `nothing` for JutulDarcy's tables (as
[`jutuldarcy_property_tables`](@ref)'s CO2 density).

Returns a new NamedTuple `(structural, residual, dissolved, trap_height,
density)`: the three masses (kg) per column, the height of each column's
trap part (m) and the CO2 density used (kg/m³). `s` and `res` are only read.
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

co2lab-demo's reachable structural capacity: for each column, the
structural capacity (`structural`, per column, e.g. from
[`trapping_capacity`](@ref)) of the trap it drains into and of every trap
downstream of that one — what CO2 injected there could fill; 0 for columns
draining out of the domain. `res` is the [`trap_analysis`](@ref) result.
Returns a new vector, one value per column.
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

# Rock compressibility: a pressure-dependent pore volume, in co2lab's form.
#
# co2lab carries it as a fluid option (makeVEFluid's `pvMult_p_ref`,
# `pvMult_fac`); JutulDarcy, like this file, keeps it with the model. The two
# forms differ slightly:
#
#   co2lab (here):   pv = pv_static · (1 + c·(p − p_ref))              exactly linear
#   JutulDarcy:      pv = pv_static · (1 + x + x²/2),  x = c·(p − p_ref)
#                    (LinearlyCompressiblePoreVolume)
#
# The x²/2 term is 5e-7 relative at 100 bar from p_ref with c = 1e-10 /Pa.
# To use JutulDarcy's form, call its function *instead of* ours, not after:
#   setup_co2lab_model(...; rock_compressibility = nothing)
#   JutulDarcy.set_rock_compressibility!(model; reference_pressure, compressibility)
# (both move the pore volume to a static parameter, under different names, so
# JutulDarcy's fails on a model ours has already changed).

"""
Pore volume `StaticFluidVolume·(1 + c·(p − p_ref))`, co2lab's `pvMultR`; set
with [`set_co2lab_rock_compressibility!`](@ref).
"""
struct LinearPoreVolume <: Jutul.ScalarVariable
    p_ref::Float64  # reference pressure (Pa)
    c::Float64      # compressibility (1/Pa)
end

# Modifies: `pv` (fluid volume).
@jutul_secondary function update_linear_pore_volume!(pv, v::LinearPoreVolume, model,
                                                    Pressure, StaticFluidVolume, ix)
    for i in ix
        pv[i] = StaticFluidVolume[i] * (1 + v.c * (Pressure[i] - v.p_ref))
    end
end

"""
    set_co2lab_rock_compressibility!(model; reference_pressure = 100e5,
                                     compressibility = 1e-10) -> model

Give the reservoir co2lab's rock compressibility: pore volume
`pv_static · (1 + c·(p − p_ref))`, exactly linear. The counterpart of
JutulDarcy's `set_rock_compressibility!` (same keywords); use one or the
other, not both. Called by [`setup_co2lab_model`](@ref).

# Arguments

- `model`: a reservoir `SimulationModel`, or a `MultiModel` (its reservoir).
- `reference_pressure`: `p_ref` (Pa); default as co2lab's `pvMult_p_ref`.
- `compressibility`: `c` (1/Pa); default as co2lab's `pvMult_fac`.

# Returns

`model` itself.

Modifies: the reservoir model's `parameters` (the static pore volume moves
from `:FluidVolume` to `:StaticFluidVolume`, once) and `secondary_variables`
(`:FluidVolume` becomes a `LinearPoreVolume`).
"""
function set_co2lab_rock_compressibility!(model; reference_pressure = 100e5,
                                          compressibility = 1e-10)
    rmodel = JutulDarcy.reservoir_model(model)
    prm = rmodel.parameters
    # the domain's pore volume becomes the static part, as JutulDarcy does
    if !haskey(prm, :StaticFluidVolume)
        prm[:StaticFluidVolume] = prm[:FluidVolume]
        delete!(prm, :FluidVolume)
    end
    rmodel.secondary_variables[:FluidVolume] =
        LinearPoreVolume(Float64(reference_pressure), Float64(compressibility))
    return model
end

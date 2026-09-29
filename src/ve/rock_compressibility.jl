# Rock compressibility: a property of the model, not of the fluid.
#
# co2lab's makeVEFluid carries it as a fluid option (pvMult_p_ref, pvMult_fac);
# JutulDarcy keeps it with the model: `set_rock_compressibility!(model;
# reference_pressure, compressibility)` (JutulDarcy/src/utils.jl:3165-3180)
# installs a pressure-dependent pore volume. This file does the same with
# co2lab's form, so the two are interchangeable (decided 2026-09-25):
#
#   co2lab (here):   pv = pv_static · (1 + c·(p − p_ref))              exactly linear
#   JutulDarcy:      pv = pv_static · (1 + x + x²/2),  x = c·(p − p_ref)
#                    (LinearlyCompressiblePoreVolume)
#
# To use JutulDarcy's native form instead, call
#   setup_co2lab_model(...; rock_compressibility = nothing)
#   JutulDarcy.set_rock_compressibility!(model; reference_pressure, compressibility)
# i.e. JutulDarcy's function *instead of* ours, not after it: both move the
# domain's pore volume into the parameter `StaticFluidVolume`, and JutulDarcy's
# looks for it under another name (`:StaticPoreVolume`, utils.jl:3173), so on
# a model where ours has already moved it, it fails looking for `:FluidVolume`.
# The two forms differ by the x²/2 term: 5e-7 relative at 100 bar from the
# reference pressure with c = 1e-10 /Pa.

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

co2lab's rock compressibility (makeVEFluid's `pvMult_p_ref`, `pvMult_fac`;
defaults as there): pore volume `pv_static · (1 + c·(p − p_ref))`, exactly
linear. The counterpart of JutulDarcy's `set_rock_compressibility!`, with the
same arguments; see the top of `rock_compressibility.jl` for exchanging the
two. `model` is a reservoir `SimulationModel` or a `MultiModel` (its
reservoir).

Modifies: the reservoir model's `parameters` (the static pore volume moves
from `:FluidVolume` to `:StaticFluidVolume`, once) and `secondary_variables`
(`:FluidVolume` becomes a `LinearPoreVolume`).

Returns `model` itself.
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

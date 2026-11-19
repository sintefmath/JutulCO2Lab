# The flux across transition faces of a VE or hybrid grid.
#
# Ordinary faces (fine–fine, continuous VE–VE) take JutulDarcy's stock
# two-point flux. A transition face (discontinuous or vertical VE–VE, VE–fine)
# treats both sides as virtual cells over the face's depth interval:
#
# - each side's phase pressures are carried hydrostatically, with the face's
#   mean phase densities, to one common depth and compared there (no separate
#   gravity term, no mixing of datums);
# - each phase takes the upwind side's mobility over its interval: a VE block's
#   relperms integrated over the interval by its vertical model, a fine cell's
#   own.
#
# The faces reach the flux through the conservation law's flux type (`VEFlux`),
# so JutulDarcy's fused two-point assembly is kept as it is.
#
# Phase pressures at depth x of cell c, with z_c its reference depth (a VE
# block's caprock, a fine cell's centroid), P its brine pressure there and Pc
# its capillary pressure: brine P + ρ_W g (x − z_c), CO2 P + Pc + ρ_G g (x − z_c),
# ρ the face's mean of the two sides' densities. For a VE block that is the
# brine continuous through the column and the CO2 at P + Δρ g h at the caprock.

"""
    VEFace{T}

One face of a VE or hybrid grid, with what a transition face's flux needs.
Built from a [`ve_domain`](@ref)'s face data; `VEFlux` holds one per face.

# Fields

- `left`: the face's left cell (block).
- `right`: the face's right cell (block).
- `category`: the face category, 0–5 as in [`ve_domain`](@ref); 3–5 are
  transition faces.
- `depth`: the common depth (m) at which a transition face compares the two
  sides' potentials; `NaN` for other faces.
- `interval`: the depth interval `(t_L, b_L, t_R, b_R)` (m) of each side's
  cells at the face, measured below that side's block top; `NaN` for
  non-transition faces.
- `zref`: the reference depths (m) of the left and right cells: a VE block's
  caprock, a fine cell's centroid.
- `is_ve`: whether the left and right cells are VE blocks.
"""
struct VEFace{T}
    left::Int                # left cell
    right::Int               # right cell
    category::Int            # face category (0–5, see `ve_domain`)
    depth::T                 # common depth of a transition (NaN otherwise)
    interval::NTuple{4, T}   # (t_L, b_L, t_R, b_R) below each block's top
    zref::NTuple{2, T}       # reference depths (caprock, or fine centroid)
    is_ve::NTuple{2, Bool}   # whether each side is a VE block
end

_is_transition(f::VEFace) = f.category in _TRANSITIONS

"""
    _ve_faces(domain) -> Union{Vector{VEFace{Float64}}, Nothing}

The faces of a [`ve_domain`](@ref) as `VEFace`s, in the order of
`domain[:neighbors]`.

Returns a new vector, or `nothing` if the domain has no transition face (the
stock two-point flux then applies throughout).
"""
function _ve_faces(domain)
    haskey(domain, :face_category) || return nothing
    cat = domain[:face_category]
    any(in(_TRANSITIONS), cat) || return nothing
    N, I, dep = domain[:neighbors], domain[:face_intervals], domain[:face_depth]
    zref, ve = domain[:reference_depth], domain[:is_ve]
    return [VEFace(N[1, f], N[2, f], cat[f], dep[f], Tuple(I[:, f]),
                   (zref[N[1, f]], zref[N[2, f]]), (ve[N[1, f]], ve[N[2, f]]))
            for f in axes(N, 2)]
end

# ------------------------------ The face's flux -------------------------------

"""
    VEFlux(faces)

The conservation law's flux type on a grid with transition faces.

JutulDarcy's fused two-point assembly passes it to `component_mass_fluxes!`
for each half-face, so transition faces get their own flux and all others
the stock one. The assembly, its speed and the sparsity pattern are
unchanged.

# Fields

- `faces`: one [`VEFace`](@ref) per face of the grid, in the domain's face
  order.
"""
struct VEFlux{T} <: Jutul.FluxType
    faces::Vector{VEFace{T}}  # every face of the grid, in the domain's order
end

"""
    _set_transition_flux!(rmodel, faces)

Install the transition flux in the reservoir model `rmodel`, by rebuilding
its mass conservation law with the flux type `VEFlux(faces)`. Returns
`rmodel`.

Modifies: `rmodel.equations[:mass_conservation]`.
"""
function _set_transition_flux!(rmodel, faces)
    law = rmodel.equations[:mass_conservation]
    N = typeof(law).parameters[4]
    rmodel.equations[:mass_conservation] =
        Jutul.ConservationLaw(law.flow_discretization, :TotalMasses, N;
                              flux = VEFlux(faces))
    return rmodel
end

"""
JutulDarcy hook: the component mass fluxes (brine; CO2 free and dissolved)
of one half-face, out of the cell `kgrad.left`. Ordinary faces use the stock
black-oil flux; transition faces compare potentials at the common depth and
take each phase's upwind mobility over the face's interval.

Returns the updated flux tuple `q` (a new value; nothing is modified).
"""
function JutulDarcy.component_mass_fluxes!(q, face, state,
                                           model::JutulDarcy.StandardBlackOilModel,
                                           ft::VEFlux, kgrad, upw)
    f = ft.faces[face]
    _is_transition(f) ||
        return JutulDarcy.component_mass_fluxes!(q, face, state, model,
                                                 Jutul.DefaultFlux(), kgrad, upw)
    ψ, λb = _transition_terms(face, state, model, f)
    # out of kgrad's left cell: the face's left-to-right flux, or its opposite
    s = kgrad.left == f.left ? 1 : -1
    ρS = reference_densities(model.system)
    rs = JutulDarcy.has_disgas(model.system) ?
         Jutul.upwind(Jutul.SPU(f.left, f.right), state.Rs, _plain(ψ[1])) : 0.0
    # brine; CO2 free and dissolved in the brine
    q = Base.setindex(q, s * ρS[1] * λb[1] * ψ[1], 1)
    q = Base.setindex(q, s * ρS[2] * (λb[2] * ψ[2] + rs * λb[1] * ψ[1]), 2)
    return q
end

# −T (Φ_R − Φ_L) per phase (brine, CO2), both at the common depth; each phase
# carried there with the face's mean density (as JutulDarcy's two-point flux),
# which keeps a compressible fluid's discrete hydrostatic equilibrium exact
function _transition_potentials(face, state, f::VEFace)
    T = state.Transmissibilities[face]
    ρ = state.PhaseMassDensities
    ρf = ntuple(ph -> (ρ[ph, f.left] + ρ[ph, f.right]) / 2, 2)
    pl = _phase_pressures_at(state, f.left, f.depth - f.zref[1], ρf)
    pr = _phase_pressures_at(state, f.right, f.depth - f.zref[2], ρf)
    return (-T * (pr[1] - pl[1]), -T * (pr[2] - pl[2]))
end

# Cell c's brine and CO2 pressures a height dz below its reference depth,
# with the phase densities ρf
function _phase_pressures_at(state, c, dz, ρf)
    p, g = state.Pressure[c], Jutul.gravity_constant
    pW = p + ρf[1] * g * dz
    pG = p + state.CapillaryPressure[1, c] + ρf[2] * g * dz
    return (pW, pG)
end

# The potentials (left to right) and each phase's upwind mobility b·kr/μ
function _transition_terms(face, state, model, f::VEFace)
    ψ = _transition_potentials(face, state, f)
    λb = ntuple(ph -> _transition_mobility(state, model, f, _plain(ψ[ph]) >= 0 ? 1 : 2,
                                           ph), 2)
    return (ψ, λb)
end

"""
    _transition_mobility(state, model, f, side, ph) -> Real

Phase `ph`'s mobility `b·kr/μ` on side `side` (1 left, 2 right) of the
transition face `f`. A fine cell gives its own; a VE block gives its relperms
integrated over the face's interval ([`_interval_relperm`](@ref)).
"""
function _transition_mobility(state, model, f::VEFace, side, ph)
    c = side == 1 ? f.left : f.right
    f.is_ve[side] || return state.SurfaceVolumeMobilities[ph, c]
    v = model.secondary_variables[:RelativePermeabilities]
    rw, rg = v.res_water, v.res_gas
    sG = state.sG[c]
    smax = _effective_max(state.sGmax[c], sG, rw, rg)
    h, hmax = state.PlumeHeights[1, c], state.PlumeHeights[2, c]
    ρ = state.PhaseMassDensities
    T, B = f.interval[2side - 1], f.interval[2side]
    kr = _interval_relperm(v.vertical, c, T, B, sG, smax, h, hmax, ρ[1, c] - ρ[2, c],
                           state.Rugosity[c], v.krw, v.krg, rw, rg)
    return state.ShrinkageFactors[ph, c] * kr[ph] / state.PhaseViscosities[ph, c]
end

# ------------------------- Relperms over an interval --------------------------

"""
    _interval_relperm(m, i, T, B, sG, smax, h, h_max, drho, dh, krw, krg, rw, rg)
        -> (krW, krG)

Column `i`'s relperms over the depths `[T, B]` (m below its top), averaged
by permeability, under the vertical model `m`. The interval `[0, H]` gives
the column's own relperms. The subscale relief `[0, dh]` takes no part in
lateral flow. Arguments may carry AD derivatives.
"""
function _interval_relperm(::SharpInterface, i, T, B, sG, smax, h, h_max, drho, dh,
                           krw, krg, rw, rg)
    return _sharp_interval_relperm(z -> z, T, B, h, h_max, dh, krw, krg)
end

function _interval_relperm(m::SharpInterfaceIntegrated, i, T, B, sG, smax, h, h_max,
                           drho, dh, krw, krg, rw, rg)
    return _sharp_interval_relperm(z -> _integral(m, i, m.K, z), T, B, h, h_max, dh,
                                   krw, krg)
end

# The sharp interface's zones within [T, B] by the permeability integral K:
# CO2 (krg) in the plume below the relief, brine 1 below h_max and krw
# between h and h_max
function _sharp_interval_relperm(K, T, B, h, h_max, dh, krw, krg)
    Tm = max(T, dh)
    z0 = zero(h) + zero(h_max)
    KTB = K(B) - K(T)
    (KTB > 0 && B > Tm) || return (z0, z0)
    clampz(z) = _min_c(_max_c(z, Tm), B)
    krG = krg * (K(clampz(h)) - K(Tm)) / KTB
    krW = (K(B) - K(clampz(h_max)) + krw * (K(clampz(h_max)) - K(clampz(h)))) / KTB
    return (krW, krG)
end

# co2lab's fringe ('P-scaled table'): krG from the table, the brine as
# krW_simple over the interval (no rugosity: refused)
function _interval_relperm(m::PScaledTable, i, T, B, sG, smax, h, h_max, drho, dh,
                           krw, krg, rw, rg)
    dg = _table_drho(m, i, drho) * Jutul.gravity_constant
    L = (B - T) * dg
    pc, pcm = h * dg, h_max * dg
    krG = (_krP(m, pc - dg * T) - _krP(m, pc - dg * B)) / L
    sf = (_SP(m, pc - dg * T) - _SP(m, pc - dg * B)) / L
    sm = (_SP(m, pcm - dg * T) - _SP(m, pcm - dg * B)) / L
    C = rg / (1 - rw)
    s = (1 - C) * sf + C * sm
    return (_max_c(1 - s - sf * rw / (1 - rw), 0.0), krG)
end

# LayeredFringe: the interval integrals below the relief, the brine with the
# tip's treatment, scaled by the interval's whole ∫ k dz
function _interval_relperm(m::LayeredFringeModel, i, T, B, sG, smax, h, h_max, drho,
                           dh, krw, krg, rw, rg)
    Tm = max(T, dh)
    z0 = zero(h) + zero(h_max) + zero(drho)
    kz = _interval_integrals(m, i, T, B, h, h_max, drho).kz
    (kz > 0 && B > Tm) || return (z0, z0)
    I = _interval_integrals(m, i, Tm, B, h, h_max, drho)
    krW = _column_brine(m.law, m, i, sG, smax, h_max, drho, Tm, B, I, rw, rg)
    return (krW / kz, I.krG / kz)
end

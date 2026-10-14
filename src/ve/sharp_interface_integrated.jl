# co2lab's 'sharp_interface_integrated' vertical model: the sharp interface in
# a column whose rock varies with depth (co2lab: addVERelpermSharpInterface
# with 'type', 'integrated'; upscaledSat2height's heterogeneous branch).
#
# - Heights by pore volume: h_max from sGmax·Vtot = (1 − rw)·V(h_max), h from
#   sG·Vtot − rg·V(h_max) = (1 − rw − rg)·V(h), with V(z) the pore volume above
#   depth z. V is piecewise linear (one piece per fine cell), so it is inverted
#   exactly (co2lab: Newton to 1e-5 m³).
# - Relperm by x-permeability: the CO2's is the plume's share of ∫ k_x dz, as in
#   co2lab. The brine's is zoned like the sharp interface's (1 below h_max, the
#   endpoint between h and h_max), where co2lab scales the whole brine column
#   by the endpoint.
# - The touched brine (dissolution) and the inventory split go by pore volume.
#
# Uniform rock gives `SharpInterface`'s results.

# --------------------------------- Building -----------------------------------

"""
    SharpInterfaceIntegrated(mesh, columns; poro, perm_x, perm_y = nothing,
                             geometry = ve_column_geometry(mesh, columns),
                             anisotropy_rtol = 1e-6) -> SharpInterfaceIntegrated

Build the integrated sharp interface for one grid, to give to
[`co2lab_ve_fluid`](@ref) as `vertical`.

# Arguments

- `mesh`, `columns`: the 3D mesh and its VE columns ([`ve_columns`](@ref)).
- `poro`, `perm_x`: porosity and x-permeability (m²) per cell of `mesh`, the
  same as given to [`ve_domain`](@ref).
- `perm_y`: only checked: the relperm uses `perm_x` only, as co2lab, and a
  warning is given if `perm_y` differs from it by more than `anisotropy_rtol`
  (relative) anywhere.
- `geometry`: the columns' geometry ([`ve_column_geometry`](@ref)), if
  already computed.

# Returns

A new `SharpInterfaceIntegrated` holding cumulative integrals per column. It
belongs to this grid: using it on a domain with a different number of
columns (e.g. in [`setup_co2lab_model`](@ref)) throws.
"""
function SharpInterfaceIntegrated(mesh, columns::VEColumns; poro, perm_x,
                                  perm_y = nothing,
                                  geometry = ve_column_geometry(mesh, columns),
                                  anisotropy_rtol = 1e-6)
    isnothing(perm_y) ||
        _warn_lateral_anisotropy(perm_x, perm_y, anisotropy_rtol;
                                 what = "the integrated relperm")
    vol = tpfv_geometry(mesh).volumes
    ncol = number_of_columns(columns)
    # one entry per cell boundary: a column of n cells has n + 1
    pos = columns.pos .+ (0:ncol)
    z, V, K = zeros(pos[end] - 1), zeros(pos[end] - 1), zeros(pos[end] - 1)
    for c in 1:ncol
        b = pos[c]
        for (j, k) in enumerate(columns.pos[c]:(columns.pos[c+1] - 1))
            cell, dz = columns.cells[k], geometry.dz[k]
            z[b+j] = z[b+j-1] + dz
            # pore volume φ·volume (co2lab: φ·volume/dz per unit depth)
            V[b+j] = V[b+j-1] + poro[cell] * vol[cell]
            K[b+j] = K[b+j-1] + perm_x[cell] * dz
        end
    end
    return SharpInterfaceIntegrated(pos, z, V, K)
end

# The model itself, after checking that the domain has its number of columns
function _vertical_model(spec::SharpInterfaceIntegrated, domain, fluid)
    n = length(domain[:column_height])
    length(spec.pos) - 1 == n ||
        error("SharpInterfaceIntegrated has $(length(spec.pos) - 1) columns, " *
              "the domain $n; build it from the domain's mesh and columns")
    return spec
end

# ------------------------------ Column integrals ------------------------------

# `Y` (V or K) of column i at depth h below the caprock (co2lab's
# integrateVertically): linear within each cell, at a boundary the slope below
function _integral(m::SharpInterfaceIntegrated, i, Y, h)
    r = m.pos[i]:(m.pos[i+1] - 1)
    k = clamp(searchsortedlast(view(m.z, r), Jutul.value(h)), 1, length(r) - 1)
    a = first(r) + k - 1
    dz = m.z[a+1] - m.z[a]
    dz > 0 || return zero(h) + Y[a+1]
    return Y[a] + (Y[a+1] - Y[a]) / dz * (h - m.z[a])
end

# The depth where `Y` (V or K) of column i reaches v, with its derivative;
# beyond the bottom extrapolated from the bottom cell (callers cap at H)
function _inverse(m::SharpInterfaceIntegrated, i, Y, v)
    r = m.pos[i]:(m.pos[i+1] - 1)
    # the cell where Y passes v; from the top cell at v = 0 (its derivative)
    k = Jutul.value(v) <= 0 ? 1 : searchsortedlast(view(Y, r), Jutul.value(v))
    a = first(r) + clamp(k, 1, length(r) - 1) - 1
    dY = Y[a+1] - Y[a]
    # a cell without pore volume: its top
    dY > 0 || return zero(v) + m.z[a]
    return m.z[a] + (v - Y[a]) * (m.z[a+1] - m.z[a]) / dY
end

_total(m::SharpInterfaceIntegrated, i, Y) = Y[m.pos[i+1] - 1]

# ------------------------- Vertical-model interface ---------------------------

# co2lab's upscaledSat2height, heterogeneous: the plume and residual zone by
# pore volume
function _column_heights(m::SharpInterfaceIntegrated, i, sG, sGmax, H, drho, rw, rg)
    Vtot = _total(m, i, m.V)
    h_max = _min_c(_inverse(m, i, m.V, sGmax * Vtot / (1 - rw)), H)
    # the residual CO2 down to h_max, the rest mobile in the plume
    v = _max_c(sG * Vtot - rg * _integral(m, i, m.V, h_max), 0.0)
    h = _min_c(_inverse(m, i, m.V, v / (1 - rw - rg)), H)
    return (h, h_max)
end

# krG_integrated (co2lab), and the brine zoned like krW_simple: relperm 1
# below h_max, krw between h and h_max, each zone weighted by its permeability
function _column_relperm(m::SharpInterfaceIntegrated, i, sG, sGmax, h, h_max, H, drho,
                         dh, krw, krg, rw, rg)
    K(x) = _integral(m, i, m.K, x)
    KH = _total(m, i, m.K)
    # the accretion layer under the caprock carries no CO2
    lost = K(min(dh, H))
    krG = _max_c(K(h) - lost, 0.0) * krg / KH
    Kmax = K(h_max)
    krW = (KH - _max_c(Kmax, lost) + _max_c(Kmax - _max_c(K(h), lost), 0.0) * krw) / KH
    return (krW, krG)
end

_fine_saturations(::SharpInterfaceIntegrated, i, zb, h, h_max, drho, rw, rg) =
    _column_saturations(zb, h, h_max, rw, rg)

# The pore volume below h_max (co2lab: the height below it, homogeneous)
function _untouched_brine(m::SharpInterfaceIntegrated, i, smax, h_max, H, drho)
    Vtot = _total(m, i, m.V)
    return (Vtot - _integral(m, i, m.V, _min_c(h_max, H))) / Vtot
end

# The inventory's split by pore volume; co2lab's by height (co2lab_compatible)
function _category_saturations(m::SharpInterfaceIntegrated, S, Smax, h, h_max, col,
                               drho, rw, rg; co2lab_compatible::Bool = false)
    co2lab_compatible && return _height_split(S, Smax, h, h_max, col, rw, rg)
    # each depth replaced by the pore volume above it: V is increasing, so the
    # height split's zones and fractions become pore-volume ones
    V(x) = [_integral(m, i, m.V, min(x[i], col.H[i])) for i in eachindex(x)]
    vcol = (H = V(col.H), zt = V(col.zt), dh = V(col.dh))
    return _height_split(S, Smax, V(h), V(h_max), vcol, rw, rg)
end

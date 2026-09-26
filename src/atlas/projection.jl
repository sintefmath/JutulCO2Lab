# UTM coordinates (transverse Mercator, GRS80 ellipsoid) by the Krüger series
# to third order in n (Karney 2011, "Transverse Mercator with an accuracy of
# a few nanometers", eqs. 35-36): millimetre accuracy within a few zones of the
# central meridian. Used to bring the coastline shipped with the Atlas (UTM
# zone 33) into the Atlas's zone 31. No datum shift: ED50 and EUREF89 differ by
# ~150 m here, below the coastline's own generalisation.

# ------------------------------- The ellipsoid --------------------------------

const _UTM_A = 6378137.0                 # GRS80 semi-major axis (m)
const _UTM_F = 1 / 298.257222101         # GRS80 flattening
const _UTM_K0 = 0.9996                   # UTM scale on the central meridian

# Third flattening n, rectifying radius A, and the series coefficients
function _utm_constants()
    n = _UTM_F / (2 - _UTM_F)
    A = _UTM_A / (1 + n) * (1 + n^2 / 4 + n^4 / 64)
    α = (n / 2 - 2n^2 / 3 + 5n^3 / 16, 13n^2 / 48 - 3n^3 / 5, 61n^3 / 240)
    β = (n / 2 - 2n^2 / 3 + 37n^3 / 96, n^2 / 48 + n^3 / 15, 17n^3 / 480)
    δ = (2n - 2n^2 / 3 - 2n^3, 7n^2 / 3 - 8n^3 / 5, 56n^3 / 15)
    return (A = A, α = α, β = β, δ = δ, e = sqrt(_UTM_F * (2 - _UTM_F)))
end

_central_meridian(zone) = 6zone - 183

# ------------------------------ Forward, inverse ------------------------------

"UTM `(x, y)` (m, northern hemisphere) of latitude, longitude (degrees)."
function _utm_forward(lat, lon, zone)
    c = _utm_constants()
    φ, λ = deg2rad(lat), deg2rad(lon - _central_meridian(zone))
    t = sinh(atanh(sin(φ)) - c.e * atanh(c.e * sin(φ)))
    ξ, η = atan(t, cos(λ)), atanh(sin(λ) / sqrt(1 + t^2))
    x = η + sum(c.α[j] * cos(2j * ξ) * sinh(2j * η) for j in 1:3)
    y = ξ + sum(c.α[j] * sin(2j * ξ) * cosh(2j * η) for j in 1:3)
    return (500000 + _UTM_K0 * c.A * x, _UTM_K0 * c.A * y)
end

"Latitude, longitude (degrees) of UTM `(x, y)` (m, northern hemisphere)."
function _utm_inverse(x, y, zone)
    c = _utm_constants()
    ξ1, η1 = y / (_UTM_K0 * c.A), (x - 500000) / (_UTM_K0 * c.A)
    ξ = ξ1 - sum(c.β[j] * sin(2j * ξ1) * cosh(2j * η1) for j in 1:3)
    η = η1 - sum(c.β[j] * cos(2j * ξ1) * sinh(2j * η1) for j in 1:3)
    χ = asin(sin(ξ) / cosh(η))
    φ = χ + sum(c.δ[j] * sin(2j * χ) for j in 1:3)
    λ = atan(sinh(η), cos(ξ))
    return (rad2deg(φ), _central_meridian(zone) + rad2deg(λ))
end

"UTM `(x, y)` of zone `from` in zone `to`; `NaN` stays `NaN`."
function _utm_rezone(x, y, from, to)
    (isnan(x) || isnan(y)) && return (NaN, NaN)
    return _utm_forward(_utm_inverse(x, y, from)..., to)
end

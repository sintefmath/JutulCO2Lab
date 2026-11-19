"""
    JutulCO2Lab.Oracle

Comparison harness for gating JutulCO2Lab against MRST-co2lab: Julia results
vs MRST results on identical input, quantity by quantity, each with a
declared tolerance, collected in a pass/fail report.

Many errors fail *silently* (the model converges and the answer is wrong),
so comparisons are numeric with explicit tolerances, not eyeballed on plots.

# Usage

    ref = read_mrst_dump("stage1_johansen.mat")
    rep = OracleReport("Johansen grid"; mrst = ref)
    check!(rep, "column height", H_julia, H_mrst; tol = STAGE1_TOL[:height])
    display(rep)
    passed(rep) || error("gate failed")

# Deliberately not done

- No guessing of cell or face orderings: MRST's and Jutul's numbering need
  not agree, and a silently permuted comparison is exactly the false pass
  this harness exists to prevent. Pass an explicit `perm` to
  [`check!`](@ref), or compare order-invariant reductions.
- No VE physics tolerances: [`STAGE1_TOL`](@ref) suits grid and geometry
  quantities, which should agree to near round-off. Flow quantities need
  looser, deliberately chosen values.
"""
module Oracle

using MAT
using Printf
using Statistics

export MrstDump, read_mrst_dump, OracleReport, Tolerance, check!, passed, STAGE1_TOL
export SourceInfo, source_provenance

include("provenance.jl")

# --------------------------------- MRST dumps ---------------------------------

"""
    MrstDump

An MRST reference dump: a `.mat` file written by the MRST/Octave dump
scripts, loaded with [`read_mrst_dump`](@ref).

Access a variable with `d["name"]`; `keys(d)` lists what the dump contains,
`haskey(d, "name")` tests for one.

# Fields

- `path`: the `.mat` file the dump was read from.
- `data`: the dump's variables by name, as MAT.jl returns them (`Float64`,
  column-major, scalars as 1×1 matrices).
"""
struct MrstDump
    path::String             # file it was read from
    data::Dict{String, Any}  # the variables, by name
end

Base.getindex(d::MrstDump, k::AbstractString) = d.data[k]
Base.haskey(d::MrstDump, k::AbstractString) = haskey(d.data, k)
Base.keys(d::MrstDump) = keys(d.data)

function Base.show(io::IO, d::MrstDump)
    print(io, "MrstDump(\"", d.path, "\", ", length(d.data), " entries)")
end

"""
    read_mrst_dump(path) -> MrstDump

Read a `.mat` file written by the MRST side of the harness.

# Arguments

- `path`: the file; v5 and the HDF5-based v7.3 format (`save(..., '-v7.3')`,
  needed for large dumps) both work.

# Returns

A new [`MrstDump`](@ref) with the file's variables as MATLAB stored them:
`Float64`, column-major, scalars as 1×1 matrices. Reshape (e.g. `vec`,
`only`) at the call site, where the expected shape is known.
"""
function read_mrst_dump(path::AbstractString)
    isfile(path) || throw(ArgumentError("no such MRST dump: $path"))
    return MrstDump(path, matread(path))
end

# --------------------------------- Tolerances ---------------------------------

"""
    Tolerance(rtol, atol)

The tolerance of one quantity's comparison in [`check!`](@ref). A value `a`
passes against the MRST reference `b` when

    |a - b| <= atol + rtol * |b|

# Fields

- `rtol`: relative tolerance, scaling the MRST value's magnitude `|b|`.
- `atol`: absolute tolerance, in the quantity's SI unit.
"""
struct Tolerance
    rtol::Float64  # relative tolerance
    atol::Float64  # absolute tolerance
end

"""
    STAGE1_TOL

Default [`Tolerance`](@ref)s for comparing grid and geometry quantities,
keyed by quantity (`:height`, `:pore_volume`, `:transmissibility`, ...).

They are tight on purpose: both sides do the same arithmetic on the same
input, so anything beyond round-off is a real discrepancy.

`atol` is in SI units, scaled to the quantity's typical magnitude (a
transmissibility of ~1e-12 m³ needs an `atol` far below a pore volume of
~1e5 m³). Pass another tolerance when a quantity's scale differs.
"""
const STAGE1_TOL = Dict{Symbol, Tolerance}(
    :height          => Tolerance(1e-12, 1e-9),    # m
    :top_depth       => Tolerance(1e-12, 1e-9),    # m
    :bottom_depth    => Tolerance(1e-12, 1e-9),    # m
    :pore_volume     => Tolerance(1e-12, 1e-6),    # m^3
    :porosity        => Tolerance(1e-12, 1e-14),   # -
    :permeability    => Tolerance(1e-12, 1e-30),   # m^2
    :net_to_gross    => Tolerance(1e-12, 1e-14),   # -
    :transmissibility=> Tolerance(1e-10, 1e-24),   # m^3 (T = k*A/d)
    :trans_mult      => Tolerance(1e-12, 1e-14),   # -
    :well_index      => Tolerance(1e-10, 1e-24),   # m^3
)

# ----------------------------- Reports and checks -----------------------------

"""
    CheckResult

The outcome of one quantity's comparison, returned by [`check!`](@ref).
`maxabs` and `maxrel` are kept even on a pass, so a gate can be seen drifting
before it breaks.

# Fields

- `name`: the quantity compared, as shown in the report.
- `n`: the number of values compared (the Julia side's length).
- `maxabs`: the largest absolute difference over the finite pairs; `NaN`
  when the comparison could not be made (length mismatch, empty).
- `maxrel`: the largest difference relative to the MRST value; `NaN` as
  for `maxabs`.
- `tol`: the [`Tolerance`](@ref) applied.
- `ok`: whether every pair was within `tol` (false for any non-finite pair,
  length mismatch or empty comparison).
- `note`: a remark shown in the report: the caller's note and/or why the
  check failed; empty if none.
"""
struct CheckResult
    name::String     # quantity compared
    n::Int           # number of values
    maxabs::Float64  # largest absolute difference
    maxrel::Float64  # largest relative difference
    tol::Tolerance   # tolerance applied
    ok::Bool         # passed
    note::String     # remark shown in the report (e.g. why it failed)
end

"""
    OracleReport(title; mrst = nothing)

A comparison report: collects [`check!`](@ref) results and displays them as
one table, with a PASS/FAIL line ([`passed`](@ref)).

# Constructor

- `title`: the report's heading.
- `mrst`: optionally the [`MrstDump`](@ref) compared against, to record the
  MRST commit it was made from.

# Fields

- `title`: the report's heading.
- `results`: one [`CheckResult`](@ref) per [`check!`](@ref), in call order.
- `sources`: the [`SourceInfo`](@ref) of Jutul and JutulDarcy, then of MRST
  if a dump was given.

The header records the commit and dirty state of the loaded Jutul and
JutulDarcy checkouts ([`source_provenance`](@ref)); construction warns if
either has uncommitted changes. See the module docstring for an example.
"""
struct OracleReport
    title::String                 # heading of the report
    results::Vector{CheckResult}  # one per check!
    sources::Vector{SourceInfo}   # commits of Jutul, JutulDarcy, MRST
end

function OracleReport(title::AbstractString; mrst = nothing)
    sources = source_provenance()
    mrst === nothing || push!(sources, _mrst_source_info(mrst))
    _warn_dirty(sources)
    return OracleReport(String(title), CheckResult[], sources)
end

"""
    passed(report) -> Bool

Whether every check in `report` passed. Returns `false` for an empty report:
a gate that ran no checks has not been executed.
"""
passed(r::OracleReport) = !isempty(r.results) && all(x -> x.ok, r.results)

"""
    check!(report, name, julia, mrst; tol, perm = nothing, note = "")

Compare `julia` against the MRST reference `mrst` elementwise and record the
result in `report`.

# Arguments

- `report`: the [`OracleReport`](@ref) to add to.
- `name`: the quantity's name in the report.
- `julia`, `mrst`: the values, scalars or arrays; arrays are flattened
  (MAT.jl returns MATLAB's 1×1 and N×1 matrices as such).
- `tol`: the [`Tolerance`](@ref).
- `perm`: maps Julia indices onto MRST indices, comparing `julia[i]` with
  `mrst[perm[i]]`. Supply it whenever the orderings differ; none is
  inferred. Default: compare in the given order.
- `note`: a remark shown in the report.

# Returns

The new [`CheckResult`](@ref). A length mismatch, an empty comparison or a
non-finite value is recorded as a failure, not thrown.

Modifies: `report.results` (appends one `CheckResult`).
"""
function check!(report::OracleReport, name::AbstractString, julia, mrst;
                tol::Tolerance,
                perm::Union{Nothing, AbstractVector{<:Integer}} = nothing,
                note::AbstractString = "")
    a = _flatten(julia)
    b = _flatten(mrst)
    # a bad perm or length is recorded as a failure, not thrown: the gate goes on
    if perm !== nothing
        if length(perm) != length(a)
            res = CheckResult(String(name), length(a), NaN, NaN, tol, false,
                              "perm length $(length(perm)) != julia length $(length(a))")
            push!(report.results, res)
            return res
        end
        b = b[perm]
    end

    if length(a) != length(b)
        res = CheckResult(String(name), length(a), NaN, NaN, tol, false,
                          "length mismatch: julia $(length(a)) vs mrst $(length(b))")
        push!(report.results, res)
        return res
    end

    if isempty(a)
        res = CheckResult(String(name), 0, NaN, NaN, tol, false, "empty comparison")
        push!(report.results, res)
        return res
    end

    maxabs, maxrel, ok, nonfinite = _compare(a, b, tol)
    if nonfinite > 0
        nf = "$nonfinite non-finite entries"
        note = isempty(note) ? nf : "$note; $nf"
    end

    res = CheckResult(String(name), length(a), maxabs, maxrel, tol, ok, String(note))
    push!(report.results, res)
    return res
end

# --------------------------- Comparison and display ---------------------------

"""
    _compare(a, b, tol) -> (maxabs, maxrel, ok, nonfinite)

Elementwise comparison of equal-length vectors `a` (Julia) and `b` (MRST).
Returns the largest absolute and relative difference over the finite pairs,
whether every pair is within `tol`, and the number of non-finite pairs (each
one a failure).
"""
function _compare(a, b, tol::Tolerance)
    maxabs = 0.0
    maxrel = 0.0
    ok = true
    nonfinite = 0
    for i in eachindex(a, b)
        ai, bi = a[i], b[i]
        if !isfinite(ai) || !isfinite(bi)
            nonfinite += 1
            ok = false
            continue
        end
        d = abs(ai - bi)
        maxabs = max(maxabs, d)
        # relative to the MRST value, guarded against division by zero
        maxrel = max(maxrel, d / max(abs(bi), eps()))
        ok &= d <= tol.atol + tol.rtol * abs(bi)
    end
    return (maxabs, maxrel, ok, nonfinite)
end

_flatten(x::Number) = [Float64(x)]
_flatten(x::AbstractArray) = vec(Float64.(x))

function Base.show(io::IO, ::MIME"text/plain", r::OracleReport)
    println(io, r.title)
    println(io, "-"^length(r.title))
    isempty(r.sources) || println(io, "  ", join(string.(r.sources), "  "))
    if isempty(r.results)
        println(io, "  (no checks run — this is a FAILURE, not a pass)")
        return
    end

    w = maximum(length(x.name) for x in r.results)
    for x in r.results
        @printf(io, "  %-3s %-*s  n=%-8d max|Δ|=%-10.3g max rel=%-10.3g",
                x.ok ? "ok" : "!!", w, x.name, x.n, x.maxabs, x.maxrel)
        if !isempty(x.note)
            print(io, "  [", x.note, "]")
        end
        println(io)
    end

    nfail = count(x -> !x.ok, r.results)
    println(io, nfail == 0 ? "  PASS ($(length(r.results)) checks)" :
                             "  FAIL ($nfail of $(length(r.results)) checks)")
end

end # module Oracle

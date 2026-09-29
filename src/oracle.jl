"""
    JutulCO2Lab.Oracle

Comparison harness for the roadmap's acceptance gates.

The pass/fail oracle is **Julia-VE vs MRST-VE on identical input** (roadmap,
"How to use this"). Every trap catalogued in feasibility §2b fails *silently* —
the model converges and the answer is wrong — so a comparison must be run
per-quantity with a declared tolerance, not eyeballed on a plot.

Two things this module deliberately does *not* do:

  - It does not guess cell/face orderings. MRST's cell and face numbering and
    Jutul's need not agree, and a silently-permuted comparison is exactly the
    kind of false pass this harness exists to prevent. Pass an explicit
    `perm` (see [`check!`](@ref)) or compare order-invariant reductions.
  - It does not define VE physics tolerances for you. [`STAGE1_TOL`](@ref) is a
    starting point for the Stage 1 *geometry* gate, where quantities should
    agree to near machine precision because both sides compute the same
    arithmetic. Flow quantities need looser, deliberately-chosen values.
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

A `.mat` file written by the MATLAB side of the harness (see
`scripts/mrst_dump_stage1.m`), loaded into a `Dict{String, Any}`.

Access with `d["name"]`; `keys(d)` lists what the dump contains.
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

Read a MATLAB `.mat` file produced by the MRST side of the harness.

MAT.jl handles both the v5 and the HDF5-based v7.3 format, so `save(..., '-v7.3')`
on the MATLAB side is fine (and is what large Johansen dumps will need).

MATLAB stores everything as `Float64` and column-major; scalars arrive as
1×1 matrices. This function leaves that as-is rather than second-guessing
each field — squeeze/`only` at the call site where the expected shape is known.

Returns a new [`MrstDump`](@ref) wrapping the file's contents.
"""
function read_mrst_dump(path::AbstractString)
    isfile(path) || throw(ArgumentError("no such MRST dump: $path"))
    return MrstDump(path, matread(path))
end

# --------------------------------- Tolerances ---------------------------------

"""
    Tolerance(rtol, atol)

Per-quantity comparison tolerance. A value passes when

    |a - b| <= atol + rtol * |b|

elementwise, with `b` the MRST reference.
"""
struct Tolerance
    rtol::Float64  # relative tolerance
    atol::Float64  # absolute tolerance
end

"""
    STAGE1_TOL

Default tolerances for the Stage 1 grid/geometry gate (roadmap Stage 1
"Gate"). These are tight on purpose: at Stage 1 both sides perform the same
arithmetic on the same input, so anything beyond round-off is a real
discrepancy, not a modelling difference.

`atol` values are in SI and are set relative to the magnitude of the quantity
(a transmissibility of ~1e-12 m³ needs an `atol` far below a pore volume of
~1e5 m³). Override per call when a quantity's scale differs from the
assumptions here.
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

Outcome of one quantity comparison. `maxabs`/`maxrel` are reported even on a
pass so that a gate can be watched drifting before it breaks.
"""
struct CheckResult
    name::String     # quantity compared
    n::Int           # number of values
    maxabs::Float64  # largest absolute difference
    maxrel::Float64  # largest relative difference
    tol::Tolerance   # tolerance applied
    ok::Bool         # passed
    note::String     # remark shown in the report
end

"""
    OracleReport(title; mrst = nothing)

Accumulates [`check!`](@ref) results and prints them as one table.

    ref = read_mrst_dump("stage1_johansen.mat")
    rep = OracleReport("Stage 1 — Johansen grid"; mrst = ref)
    check!(rep, "column height", H_julia, H_mrst; tol = STAGE1_TOL[:height])
    display(rep)
    passed(rep) || error("gate failed")

On construction the report records the commit and dirty state of the loaded
Jutul and JutulDarcy checkouts (see [`source_provenance`](@ref)), and warns if
either has uncommitted changes. Pass the MRST dump as `mrst` to record the MRST
commit it was generated from as well. All of it is printed in the header.
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

True when every check in `report` passed. An empty report is *not* a pass —
a gate that ran no checks has not been executed.
"""
passed(r::OracleReport) = !isempty(r.results) && all(x -> x.ok, r.results)

"""
    check!(report, name, julia, mrst; tol, perm = nothing, note = "")

Compare `julia` against the MRST reference `mrst` elementwise and record the
result in `report`. Returns the [`CheckResult`](@ref).

Modifies: `report.results` (appends one `CheckResult`). `julia` and `mrst` are
only read.

`perm` maps Julia indices onto MRST indices: the comparison is
`julia[i]` vs `mrst[perm[i]]`. Supply it whenever the two orderings differ —
there is no attempt to infer one. Without it, the arrays are compared in
their given order and a length mismatch is a failure, not an error.

Scalars and arrays are both accepted; MATLAB's 1×1 and N×1 matrices are
flattened first, since that is how MAT.jl returns them.
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

Elementwise comparison of equal-length vectors `a` (Julia) and `b` (MRST):
largest absolute and relative difference over the finite pairs, whether every
pair is within `tol`, and the number of non-finite pairs (each one a failure).
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

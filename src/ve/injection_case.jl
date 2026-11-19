# co2lab-demo's exploreSimulation scenario: CO2 injected from several wells
# into a VE formation, then left to migrate, each boundary side closed,
# semi-open or open; and a run that hands over each finished report step
# (for a live display).

# MRST's year (365.2425 days), in which exploreSimulation gives its times
const MRST_YEAR = 365.2425 * 86400

"""
    setup_injection_case(a, wells; fluid, boundary = :open, inj_time = 50 years,
                         inj_steps = 10, mig_time = 3000 years, mig_steps = 30,
                         temperature = nothing, rugosity = 0.0, rho_w = nothing,
                         well_radius = 0.3, outside_distance = 100e3,
                         pore_compressibility = 1e-10, dissolved = :carry,
                         well_volume = 1.0, co2lab_compatible = false) -> NamedTuple

co2lab-demo's `exploreSimulation` launch, without its GUI: CO2 injected into
a VE formation from several wells, then left to migrate. Run it with
[`simulate_injection_case`](@ref) or `simulate_reservoir`.

# Arguments

- `a`: a formation as [`atlas_ve_domain`](@ref) returns it (`mesh`,
  `columns`, `ve`).
- `wells`: a vector of `(column, rate)` pairs, `rate` the CO2 mass rate
  (kg/s).
- `fluid`: the fluid, from [`co2lab_ve_fluid`](@ref).
- `boundary`: `:closed`, `:open` or `:semi_open`, one for all sides or one
  per side of [`ve_boundary_sides`](@ref)'s `sides` (below).
- `inj_time`, `inj_steps`: injection time (s; default 50 years) in equal
  report steps.
- `mig_time`, `mig_steps`: migration time (s; default 3000 years) with the
  wells shut, in equal report steps.
- `temperature` (K; a number, one per column or a function of depth, as
  [`ve_domain`](@ref)) and `rugosity` (m): replace the domain's.
- `rho_w`: brine density (kg/m³) for co2lab's hydrostatic pressure
  `rho_w·g·z`; default `nothing`, the fluid's brine at rest
  ([`ve_hydrostatic_profile`](@ref), at the formation's temperatures).
- `well_radius`: radius (m) of co2lab's VE wells ([`setup_ve_well`](@ref)).
- `outside_distance`: length (m) of the imagined formation behind a
  semi-open side.
- `pore_compressibility`: (1/Pa) about the mean initial pressure (co2lab's
  `pvMult`).
- `dissolved`: how brine crossing an open side books its dissolved CO2
  ([`ve_boundary_conditions`](@ref); `:co2lab` for comparisons).
- `well_volume`: wellbore volume as a multiple of JutulDarcy's default.
- `co2lab_compatible`: co2lab's face density, for comparisons with co2lab
  (as [`setup_co2lab_model`](@ref)). Default `false`.
  co2lab's wells have none; `1e-3` comes close, for comparisons.

# Returns

A new NamedTuple:

- `case`: the `JutulCase`;
- `model`: its model;
- `ve`: the formation `ve`, its domain copied and modified;
- `sides`, `side_types`: the boundary sides and their types (`Symbol`s);
- `bc`: the boundary conditions;
- `wells`: the well columns;
- `dt`: the report steps (s);
- `pressure0`: the initial pressure (Pa).

# Boundary

- Sealed fault sides (`across_fault`) are closed whatever is given.
- Open sides are held at hydrostatic brine pressure at the side's depth,
  with brine flowing in.
- A semi-open side is an open one behind an imagined stretch of formation
  `outside_distance` long: its transmissibility is scaled by
  `d/(d + outside_distance)`, `d` the map distance from the centroid of the
  column's top face to the side (co2lab's rule).

# Notes

- The initial state is hydrostatic brine pressure at the caprock, no CO2.
- A well whose pressure crosses CO2's boiling curve (a shallow, cold
  caprock) converges slowly: JutulDarcy's well injects at the density of its
  own pressure, which jumps there (co2lab converts at the cell's).
"""
function setup_injection_case(a, wells; fluid, boundary = :open,
                              inj_time = 50MRST_YEAR, inj_steps::Integer = 10,
                              mig_time = 3000MRST_YEAR, mig_steps::Integer = 30,
                              temperature = nothing, rugosity = 0.0, rho_w = nothing,
                              well_radius = 0.3, outside_distance = 100e3,
                              pore_compressibility = 1e-10, dissolved::Symbol = :carry,
                              well_volume = 1.0, co2lab_compatible::Bool = false)
    ve = _injection_formation(a.ve, temperature, rugosity)
    g = Jutul.gravity_constant
    sides, _ = ve_boundary_sides(a.mesh, a.columns, ve.domain)
    types = _side_types(sides, boundary)
    # hydrostatic brine: the fluid's at rest, or co2lab's constant rho_w
    p_of_z = isnothing(rho_w) ? ve_hydrostatic_profile(fluid, ve; g) : z -> rho_w * g * z
    bc = _injection_bc(a, ve, sides, types, p_of_z, outside_distance, dissolved)

    cols = [Int(w[1]) for w in wells]
    vewells = [setup_ve_well(ve, c; name = Symbol("I$k"), radius = well_radius,
                             volume_multiplier = well_volume)
               for (k, c) in enumerate(cols)]
    p0 = p_of_z.(ve.geometry.z_top)
    # the pore volume's reference pressure: the mean initial pressure
    rock = (reference_pressure = sum(p0) / length(p0),
            compressibility = pore_compressibility)
    model = setup_co2lab_model(ve.domain; wells = vewells, fluid,
                               rock_compressibility = rock, co2lab_compatible)
    # injection: each well at its mass rate, in the fluid's reference volumes
    rate(w) = TotalRateTarget(w[2] / fluid.rho_ref[2])
    inj = Dict(Symbol("I$k") => co2_injector_control(rate(w), fluid)
               for (k, w) in enumerate(wells))
    shut = Dict(Symbol("I$k") => DisabledControl() for k in eachindex(wells))
    f_inj = setup_reservoir_forces(model; control = inj, bc)
    f_shut = setup_reservoir_forces(model; control = shut, bc)
    dt = [fill(inj_time / inj_steps, inj_steps); fill(mig_time / mig_steps, mig_steps)]
    forces = [k <= inj_steps ? f_inj : f_shut for k in eachindex(dt)]
    state0 = setup_co2lab_state(model; pressure = p0, forces = f_inj)
    case = JutulCase(model, dt, forces; state0, parameters = setup_parameters(model))
    return (case = case, model = model, ve = ve, sides = sides, side_types = types,
            bc = bc, wells = cols, dt = dt, pressure0 = p0)
end

# The formation `ve` with its domain copied, the temperature and rugosity set
function _injection_formation(ve, temperature, rugosity)
    domain = deepcopy(ve.domain)
    n = length(ve.geometry.z_top)
    if !isnothing(temperature)
        T = temperature isa Function ? temperature.(ve.geometry.z_top) :
            temperature isa Number ? fill(Float64(temperature), n) : temperature
        domain[:temperature, Cells()] = collect(Float64, T)
    end
    domain[:rugosity, Cells()] = rugosity isa Number ? fill(Float64(rugosity), n) :
                                 collect(Float64, rugosity)
    return (; ve..., domain = domain)
end

# Each side's type: `boundary` for all or per side; sealed fault sides closed
function _side_types(sides, boundary)
    types = boundary isa Symbol ? fill(boundary, length(sides)) :
            collect(Symbol, boundary)
    length(types) == length(sides) ||
        throw(ArgumentError("one boundary type per side needed"))
    all(in((:closed, :open, :semi_open)), types) ||
        throw(ArgumentError("boundary types: :closed, :open or :semi_open"))
    types[sides.across_fault] .= :closed
    return types
end

# The conditions of the open and semi-open sides: hydrostatic pressure at the
# side's depth (`p_of_z`), a semi-open side's transmissibility scaled
# by d/(d + L), d the map distance from the centroid of its column's top face
# to the side
function _injection_bc(a, ve, sides, types, p_of_z, L, dissolved)
    ix = findall(!=(:closed), types)
    cx, cy = _top_face_centroids(a.mesh, a.columns)
    trans = map(ix) do k
        c, (ex, ey) = sides.column[k], sides.edge_midpoint[k]
        d = hypot(cx[c] - ex, cy[c] - ey)
        sides.half_trans[k] * (types[k] == :semi_open ? d / (d + L) : 1.0)
    end
    return ve_boundary_conditions(ve, sides, ix; pressure = p_of_z, trans,
                                  dissolved)
end

# Map-view centroid of each column's top face: the trap surface's (co2lab's
# Gt.cells.centroids) for the formation's top, the raw top face's for columns
# below it (stacked or hybrid grids)
function _top_face_centroids(mesh, columns)
    n = number_of_columns(columns)
    cx, cy = fill(NaN, n), fill(NaN, n)
    s = trap_surface(mesh, columns)
    cx[s.blocks], cy[s.blocks] = s.cx, s.cy
    rest = setdiff(1:n, s.blocks)
    isempty(rest) && return (cx, cy)
    t = ve_top_surface(mesh, columns; blocks = rest)
    _, x, y = _polygon_geometry(t.x, t.y, t.polygons)
    cx[rest], cy[rest] = x, y
    return (cx, cy)
end

"""
    simulate_injection_case(c; on_report = nothing, kwarg...) -> result

Simulate an injection case with co2lab's simulator settings, handing over
each report step as soon as it is finished (e.g. for a live display). The
results are those of one `simulate_reservoir` run.

# Arguments

- `c`: the case from [`setup_injection_case`](@ref).
- `on_report`: `nothing`, or a function `on_report(k, state)` called with
  the report step's number and its output state (a `Dict` per model, as
  `simulate_reservoir`'s states). Returning `false` stops the run.
- `kwarg...`: passed to [`setup_co2lab_simulator`](@ref) (e.g. `info_level`).

# Returns

`simulate_reservoir`'s result (`ReservoirSimResult`), or `nothing` if
`on_report` stopped the run.
"""
function simulate_injection_case(c; on_report = nothing, kwarg...)
    sim, cfg = setup_co2lab_simulator(c.case; kwarg...)
    if !isnothing(on_report)
        ends = cumsum(c.dt)
        progress = (t = Ref(0.0), k = Ref(1))
        cfg[:post_ministep_hook] = (done, report, sim, dt, forces, max_iter, cfg) ->
            _report_hook!(progress, ends, on_report, done, report, sim, dt)
    end
    try
        return simulate_reservoir(c.case; simulator = sim, config = cfg)
    catch err
        err isa _StopRun && return nothing
        rethrow()
    end
end

# Thrown by the report hook to stop a run (on_report returned false)
struct _StopRun <: Exception end

# The hook after a ministep: a converged one advances the time; at the end of
# report step k, calls on_report(k, state) (false: stop the run). Returns
# (done, report) for Jutul.
# Modifies: `progress` (time and next report step).
function _report_hook!(progress, ends, on_report, done, report, sim, dt)
    done || return (done, report)
    progress.t[] += dt
    k = progress.k[]
    if k <= length(ends) && progress.t[] >= ends[k] * (1 - 1e-12)
        progress.k[] = k + 1
        go = on_report(k, Jutul.get_output_state(sim))
        go === false && throw(_StopRun())
    end
    return (done, report)
end

"""
    ve_boundary_loops(s, sides) -> Vector{Vector{Int}}

The outer sides of a formation in order around each boundary loop (the
outer edge and each hole), so that a stretch of boundary is a range of a
loop; co2lab's `find_boundary_loops` (exploreSimulation).

# Arguments

- `s`: the surface ([`trap_surface`](@ref) or [`ve_top_surface`](@ref)) whose
  `outline` edges the sides lie on; sides are matched to edges by their
  midpoints.
- `sides`: the boundary sides ([`ve_boundary_sides`](@ref)).

# Returns

A new vector of loops, each a vector of indices into `sides`, consecutive
sides sharing a corner. Sides that match no outline edge are left out.
"""
function ve_boundary_loops(s, sides)
    # each side's outline edge, by the edge's midpoint
    key(x, y) = (round(x; digits = 3), round(y; digits = 3))
    edge_at = Dict(key((s.x[a] + s.x[b]) / 2, (s.y[a] + s.y[b]) / 2) => (a, b)
                   for (a, b) in s.outline)
    ends = Dict{Int, Tuple{Int, Int}}()
    for k in 1:length(sides)
        e = get(edge_at, key(sides.edge_midpoint[k]...), nothing)
        isnothing(e) || (ends[k] = e)
    end
    # the sides at each corner
    at_node = Dict{Int, Vector{Int}}()
    for (k, (a, b)) in ends, n in (a, b)
        push!(get!(at_node, n, Int[]), k)
    end
    loops, done = Vector{Int}[], Set{Int}()
    for start in sort!(collect(keys(ends)))
        start in done && continue
        loop, k, n = Int[], start, ends[start][2]
        # from side to side over shared corners, until back at the start
        while !(k in done)
            push!(loop, k); push!(done, k)
            next = findfirst(j -> !(j in done), at_node[n])
            isnothing(next) && break
            k = at_node[n][next]
            n = ends[k][1] == n ? ends[k][2] : ends[k][1]
        end
        push!(loops, loop)
    end
    return loops
end

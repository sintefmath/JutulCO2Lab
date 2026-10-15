# co2lab-demo's exploreSimulation: set up a CO2
# injection into a formation of the CO2 Storage Atlas by clicking — wells,
# and which stretches of the boundary are closed, semi-open or open — then run
# it (`setup_injection_case`, `simulate_injection_case`) and watch the plume
# grow and migrate, step by step, with the trapping inventory at the end.
#
# The setup is the interactive trap viewer (`interactive_traps`) with this
# app's widgets on its panel; the results have a figure of their own beside
# it (in a browser the page holds both; GLMakie opens a window for each).
#
#   julia --project=examples examples/demos/explore_simulation.jl
#   JUTULCO2LAB_BACKEND=WGLMakie julia --project=examples \
#       examples/demos/explore_simulation.jl   (then open http://localhost:8080)
#
# Clicks on the caprock, by the mode chosen: "Place wells" — a left click adds
# a well there, or selects the well nearby; the rate slider sets the selected
# well's rate; a right click removes the nearest well. "Edit boundary" — a left
# click near the boundary marks the start of a stretch, a second its end: the
# stretch between gets the boundary type chosen (red closed, orange
# semi-open, blue open; "Set all": the whole boundary; a right click cancels).
# The panel shows what the next click does. The fluid
# is JutulDarcy's CO2–brine (tables in pressure and temperature, the
# temperature from 7 °C at the sea floor, 100 m, and 35.6 °C/km); residual
# saturations 0.11 and 0.21; dissolution at co2lab's rate, at most 53 kg CO2
# per m³ brine; the capillary fringe co2lab's 'P-scaled table'; subscale
# trapping a rugosity of the caprock (m), not with the fringe (either toggle
# switches the other off). A well in a shallow, cold column
# (where the injected CO2 crosses its boiling curve) makes the run slow.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Printf

# Work started by a widget or the mouse runs as a task: in the browser a
# callback holds up the messages the plots need (errors in a task are shown,
# not lost)
background(f) = @async try
    f()
catch err
    @error "explore simulation" exception = (err, catch_backtrace())
end

const YEAR = 365.2425 * 86400
const BOUNDARY_TYPES = [:closed, :semi_open, :open]
const BOUNDARY_COLOURS = Dict(:closed => RGBf(0.85, 0, 0),
                              :semi_open => RGBf(1, 0.6, 0),
                              :open => RGBf(0.1, 0.5, 1))

# ------------------------------- The formations -------------------------------

# Formations read, by (formation, coarsening), shared by all pages: the
# formation, its trap surface, traps, boundary sides and loops
const FORMATIONS = Dict{Tuple{String, Int}, Any}()
const FORMATIONS_LOCK = ReentrantLock()

# A formation, from the cache if read before; without Atlas rock values,
# porosity 0.2 and 100 mD (co2lab's defaults)
function formation(name, c)
    lock(FORMATIONS_LOCK) do
        get!(FORMATIONS, (name, c)) do
            a = try
                atlas_ve_domain(name; coarsening = c)
            catch err
                err isa ArgumentError || rethrow()
                atlas_ve_domain(name; coarsening = c, poro = 0.2, perm = 1e-13)
            end
            s = trap_surface(a.mesh, a.columns)
            sides, _ = ve_boundary_sides(a.mesh, a.columns, a.ve.domain)
            (a = a, s = s, traps = trap_analysis(s; method = :node), sides = sides,
             loops = ve_boundary_loops(s, sides))
        end
    end
end

# ---------------------------------- The setup ---------------------------------

# This app's widgets in layout g (the top of the viewer's panel). Returns them.
# (Makie.Button, Makie.Slider: WGLMakie also exports Bonito's, by those names.)
function setup_widgets(g, name, c)
    menu(pos, options, default; w = 150) = Menu(pos; options, default,
                                                searchable = false, width = w)
    Label(g[1, 1], "Formation"; halign = :left)
    f = menu(g[2, 1], atlas_grid_formations(), name; w = 190)
    Label(g[1, 2], "Coarsening"; halign = :left)
    # (an integer default is a position in the options: here c)
    k = menu(g[2, 2], 1:6, c; w = 100)
    Label(g[3, 1], "Clicks"; halign = :left)
    mode = menu(g[4, 1], ["Place wells", "Edit boundary"], "Place wells"; w = 190)
    Label(g[3, 2], "Boundary"; halign = :left)
    btype = menu(g[4, 2], ["closed", "semi-open", "open"], "open"; w = 100)
    # what a click does, in the mode chosen
    help = Label(g[5, 1:2], ""; halign = :left, justification = :left,
                 tellwidth = false)
    tools = GridLayout(g[6, 1:2])
    set_all = Makie.Button(tools[1, 1]; label = "Set all")
    clear = Makie.Button(tools[1, 2]; label = "Clear wells")
    sl = SliderGrid(g[7, 1:2],
                    (label = "Rate (Mt/yr)", range = 0.1:0.1:10, startvalue = 1.0),
                    (label = "Rugosity (m)", range = 0:0.1:2, startvalue = 0.5);
                    tellwidth = false)
    return (formation = f, coarsening = k, mode = mode, btype = btype, help = help,
            set_all = set_all, clear = clear, rate = sl.sliders[1],
            rugosity = sl.sliders[2], run = run_widgets(g[8, 1:2]))
end

# The run's widgets: times, options, launch and stop, the status
function run_widgets(pos)
    g = GridLayout(pos)
    box(r, c, label, value) = (Label(g[r, c], label; halign = :left);
                               Textbox(g[r, c + 1]; stored_string = value, width = 60,
                                       validator = Float64))
    inj, inj_n = box(1, 1, "Injection (yr)", "50"), box(1, 3, "steps", "10")
    mig, mig_n = box(2, 1, "Migration (yr)", "3000"), box(2, 3, "steps", "30")
    opts = GridLayout(g[3, 1:4])
    toggle(c, label) = (Label(opts[1, 2c - 1], label; halign = :left);
                        Toggle(opts[1, 2c]; active = false))
    diss, fringe, subscale = toggle(1, "Dissolution"), toggle(2, "Fringe"),
                             toggle(3, "Subscale")
    buttons = GridLayout(g[4, 1:4])
    launch = Makie.Button(buttons[1, 1]; label = "Launch simulation")
    stop = Makie.Button(buttons[1, 2]; label = "Stop")
    status = Label(g[5, 1:4], ""; halign = :left, justification = :left,
                   tellwidth = false)
    return (inj = inj, inj_n = inj_n, mig = mig, mig_n = mig_n, dissolution = diss,
            fringe = fringe, subscale = subscale, launch = launch, stop = stop,
            status = status)
end

# The setup's state: the formation shown, the wells (column, rate in Mt/yr),
# the well selected (0: none), the boundary types per side, the side of the
# first click of a boundary stretch (0: none), the overlay plots
new_setup(f) = (f = Ref{Any}(f), wells = Tuple{Int, Float64}[], selected = Ref(0),
                types = Ref(fill(:open, length(f.sides))), pending = Ref(0),
                plots = Ref{Any}(nothing))

# The wells and boundary types on the caprock: made with each drawing
# (on_draw), then given new data (making plots in a browser is slow)
function overlay_plots!(v, st)
    s, δ = st.f[].s, 3e-3 * -(reverse(extrema(st.f[].s.z))...)
    seg = side_segments(s, st.f[].sides, δ)
    # thick, so that a stretch's type shows at the whole formation's scale
    bnd = lines!(v.ax, Point3f[]; color = RGBf[], linewidth = 8)
    wells = linesegments!(v.ax, Point3f[]; color = RGBf[], linewidth = 3)
    # (outlined in white: readable on any colour of the caprock)
    labels = text!(v.ax, Point3f[]; text = String[], color = RGBf[], fontsize = 15,
                   strokecolor = :white, strokewidth = 2)
    # the start of a boundary stretch, once clicked
    start = scatter!(v.ax, Point3f[]; color = :magenta, markersize = 20,
                     strokecolor = :white, strokewidth = 2)
    st.plots[] = (bnd = bnd, wells = wells, labels = labels, start = start,
                  seg = seg, δ = δ)
    update_overlays!(v, st)
end

# Each boundary side as a line along its outline edge (matched by the edge's
# midpoint, as ve_boundary_loops), just above the caprock: two points per
# side. Sides without an outline edge: a line across, a column wide.
function side_segments(s, sides, δ)
    key(x, y) = (round(x; digits = 3), round(y; digits = 3))
    edge_at = Dict(key((s.x[a] + s.x[b]) / 2, (s.y[a] + s.y[b]) / 2) => (a, b)
                   for (a, b) in s.outline)
    pts = Point3f[]
    for (k, (x, y)) in enumerate(sides.edge_midpoint)
        e = get(edge_at, key(x, y), nothing)
        if isnothing(e)
            z, n = sides.edge_depth[k] - δ, sides.normal[k]
            l = sqrt(s.area[sides.column[k]]) / 2
            push!(pts, Point3f(x - n[2] * l, y + n[1] * l, z),
                  Point3f(x + n[2] * l, y - n[1] * l, z))
        else
            push!(pts, (Point3f(s.x[i], s.y[i], s.z[i] - δ) for i in e)...)
        end
    end
    return pts
end

# The boundary as polylines around its loops (`seg` from side_segments), one
# per stretch of a type, NaN between them; and a colour per point
function boundary_lines(loops, seg, types)
    pts, cols = Point3f[], RGBf[]
    put!(q, c) = (push!(pts, q); push!(cols, c))
    gap = Point3f(NaN, NaN, NaN)
    for loop in loops, (i, k) in enumerate(loop)
        a, b = seg[2k - 1], seg[2k]
        c = BOUNDARY_COLOURS[types[k]]
        # each side from the corner it shares with the side before
        prev = i == 1 ? gap : pts[end]
        a == prev || b != prev || ((a, b) = (b, a))
        # a new polyline at a loop's start, a change of type or a break
        if i == 1 || a != prev || types[k] != types[loop[i - 1]]
            isempty(pts) || put!(gap, c)
            put!(a, c)
        end
        put!(b, c)
    end
    return pts, cols
end

# The overlays' data: the boundary, the start of a stretch, the wells as lines
# above their columns with their numbers (the selected one magenta); the help
# and the wells' summary
function update_overlays!(v, st)
    p, s, f = st.plots[], st.f[].s, st.f[]
    isnothing(p) && return
    pts, cols = boundary_lines(f.loops, p.seg, st.types[])
    Makie.update!(p.bnd; arg1 = pts, color = cols)
    k = st.pending[]
    Makie.update!(p.start; arg1 = k == 0 ? Point3f[] :
                  [Point3f(f.sides.edge_midpoint[k]..., f.sides.edge_depth[k] - p.δ)])
    zmin, zmax = extrema(s.z)
    top(c) = Point3f(s.cx[c], s.cy[c], zmin - 0.1 * (zmax - zmin))
    pts = Point3f[q for (c, _) in st.wells
                  for q in (Point3f(s.cx[c], s.cy[c], s.zc[c]), top(c))]
    colour(i) = i == st.selected[] ? RGBf(0.9, 0, 0.9) : RGBf(0, 0, 0)
    Makie.update!(p.wells; arg1 = pts,
                  color = [colour(i) for i in eachindex(st.wells) for _ in 1:2])
    Makie.update!(p.labels; arg1 = Point3f[top(c) for (c, _) in st.wells],
                  text = ["I$i" for i in eachindex(st.wells)],
                  color = colour.(eachindex(st.wells)))
    v.widgets.help.text[] = help_text(v, st)
    v.widgets.run.status.text[] = wells_text(st)
end

# What the next click does, in the mode chosen
function help_text(v, st)
    if v.widgets.mode.selection[] == "Place wells"
        return "Left click: add a well, or select one nearby (the rate\n" *
               "slider then sets its rate). Right click: remove the nearest."
    end
    st.pending[] == 0 &&
        return "Left click near the boundary: the start of a stretch\n" *
               "(red closed, orange semi-open, blue open)."
    return "Left click: the end of the stretch (magenta: its start);\n" *
           "it becomes $(v.widgets.btype.selection[]). Right click: cancel."
end

# The wells in two lines (the panel's height is fixed): their number and
# total rate, the selected one
function wells_text(st)
    isempty(st.wells) && return "No wells: click the caprock (Place wells)."
    total = sum(last, st.wells)
    head = @sprintf("%d well%s, %.1f Mt/yr in all", length(st.wells),
                    length(st.wells) == 1 ? "" : "s", total)
    k = st.selected[]
    k == 0 && return head * "\nClick a well to select it."
    return head * @sprintf("\nSelected: I%d, column %d, %.1f Mt/yr", k,
                           st.wells[k]...)
end

# A click on column `col` (0: none), by the mode chosen
function on_click!(v, st, col, button)
    col == 0 && return
    w = v.widgets
    if w.mode.selection[] == "Place wells"
        place_well!(st, col, button, w.rate)
    elseif button == Mouse.right
        st.pending[] = 0
    else
        edit_boundary!(st, col, BOUNDARY_TYPES[w.btype.i_selected[]])
    end
    update_overlays!(v, st)
end

# A well click on column `col`: the right button removes the nearest well; the
# left selects a well within two column widths, else adds one (at most ten) at
# the rate slider's value, selected. A selected well's rate shown on the slider.
function place_well!(st, col, button, rate)
    s = st.f[].s
    d(c) = hypot(s.cx[c] - s.cx[col], s.cy[c] - s.cy[col])
    near = isempty(st.wells) ? 0 : argmin([d(c) for (c, _) in st.wells])
    if button == Mouse.right
        near == 0 && return
        deleteat!(st.wells, near)
        st.selected[] = 0
    elseif near > 0 && d(st.wells[near][1]) <= 2 * sqrt(s.area[col])
        st.selected[] = near
        set_close_to!(rate, st.wells[near][2])
    elseif length(st.wells) < 10
        push!(st.wells, (col, rate.value[]))
        st.selected[] = length(st.wells)
    end
end

# The boundary stretch between two clicks (the sides nearest them, along the
# loop the first lies on) set to `type`
function edit_boundary!(st, col, type)
    f = st.f[]
    s, sides = f.s, f.sides
    d(k) = hypot((sides.edge_midpoint[k] .- (s.cx[col], s.cy[col]))...)
    near = argmin(d.(1:length(sides)))
    if st.pending[] == 0
        st.pending[] = near
        return
    end
    first_side, st.pending[] = st.pending[], 0
    loop = only(filter(l -> first_side in l, f.loops))
    near in loop || return
    i, j = findfirst(==(first_side), loop), findfirst(==(near), loop)
    # the shorter way round the loop
    n = length(loop)
    span = mod(j - i, n) <= n ÷ 2 ? (i:i + mod(j - i, n)) : (j:j + mod(i - j, n))
    st.types[][loop[mod1.(span, n)]] .= type
end

# ---------------------------------- The results -------------------------------

# The variables shown: name, value per column of a state (`p0`: the initial
# pressure), colour map
const VARIABLES = [
    ("CO2 saturation", (r, p0) -> r[:sG], :viridis),
    ("pressure (MPa)", (r, p0) -> r[:Pressure] ./ 1e6, :viridis),
    ("overpressure (MPa)", (r, p0) -> (r[:Pressure] .- p0) ./ 1e6, :plasma),
    ("plume height (m)", (r, p0) -> r[:PlumeHeights][2, :], :viridis),
    ("max CO2 saturation", (r, p0) -> r[:sGmax], :viridis),
    ("dissolved CO2 (rs)", (r, p0) -> get(r, :Rs, zero(r[:sG])), :viridis)]

# The results figure: the caprock in 3D, the variable, the time, the status
# and the trapping inventory; the last run and whether it is going on.
# Returns a NamedTuple.
function results_figure(size)
    fig = Figure(; size)
    # :free zooms the whole scene (as the setup viewer), so that close-ups are
    # not clipped at the axis box
    ax = Axis3(fig[1, 1]; zreversed = true, viewmode = :free, clip = false,
               xlabel = "x (m)", ylabel = "y (m)",
               zlabel = "depth (m)", xticks = LinearTicks(3), yticks = LinearTicks(3),
               title = "Launch a simulation")
    # the inventory under the caprock; the controls and colour bar beside both
    inv = Axis(fig[2, 1]; height = 260, title = "Trapping inventory",
               xlabel = "years", ylabel = "Mt")
    side = GridLayout(fig[1:2, 2]; tellheight = false, valign = :top)
    colsize!(fig.layout, 2, Fixed(240))
    var = Menu(side[1, 1]; options = first.(VARIABLES), searchable = false, width = 220)
    time = Makie.Slider(side[2, 1]; range = 0:1, startvalue = 0, width = 220)
    when = Label(side[3, 1], ""; halign = :left)
    status = Label(side[4, 1], ""; halign = :left, justification = :left)
    cbar = Colorbar(side[5, 1]; limits = (0.0, 1.0), height = 300, halign = :left)
    return (fig = fig, ax = ax, side = side, var = var, time = time, when = when,
            status = status, inv = inv, cbar = cbar, surface = Ref{Any}(nothing),
            states = Any[], run = Ref{Any}(nothing), running = Ref(false))
end

# the Makie extension's surface triangulation (the column of each mesh point),
# to recolour the caprock without making it anew
const EXT = Base.get_extension(JutulCO2Lab, :JutulCO2LabMakieExt)

# The results' caprock for surface `s`, coloured later (show_state!)
function results_surface!(r, s)
    empty!(r.ax)
    plot_ve_surface!(r.ax, s; color = zeros(length(s.zc)), strokewidth = 0,
                     shading = true)
    r.surface[] = (plot = r.ax.scene.plots[end],
                   vcols = EXT._surface_triangles(s, 0.0)[3])
    span(x) = -(reverse(extrema(x))...)
    m = max(span(s.x), span(s.y))
    r.ax.aspect[] = (span(s.x) / m, span(s.y) / m, 10 * max(span(s.z), 1.0) / m)
    plot_ve_outline!(r.ax, s)
end

# The run's report step k (time slider) in the variable chosen
function show_state!(r)
    run, k = r.run[], r.time.value[]
    (isnothing(run) || k < 1 || k > length(r.states)) && return
    name, value, cmap = VARIABLES[r.var.i_selected[]]
    x = Float64.(value(r.states[k], run.case.pressure0))
    lims = extrema(x)
    lims = lims[2] > lims[1] ? lims : (lims[1], lims[1] + 1)
    sp = r.surface[]
    sp.plot.colormap[] = cmap
    sp.plot.colorrange[] = lims
    sp.plot.color[] = x[sp.vcols]
    r.cbar.limits[] = lims; r.cbar.colormap[] = cmap; r.cbar.label[] = name
    r.when.text[] = @sprintf("year %.0f (step %d of %d)", cumsum(run.case.dt)[k] / YEAR,
                             k, length(run.case.dt))
end

# A run from the setup: the case, then the simulation in the background,
# each report step shown as it comes, the trapping inventory at the end
function launch!(v, st, r)
    w = v.widgets.run
    r.running[] && (r.status.text[] = "A run is going on."; return)
    isempty(st.wells) && (r.status.text[] = "Place at least one well."; return)
    num(t) = parse(Float64, t.stored_string[])
    f = st.f[]
    fluid = co2lab_ve_fluid(; vertical = w.fringe.active[] ?
                                         CapillaryFringeTable(invPc3D = (0.1, 0.5)) :
                                         SharpInterface(),
                            residual = (0.11, 0.21), dissolution = w.dissolution.active[],
                            dis_rate = 3e-11, dis_max_mass = 0.053)
    # the caprock's temperature: 7 °C at the sea floor (100 m), 35.6 °C/km
    T(z) = 273.15 + 7 + (z - 100) * 0.0356
    wells = [(col, rate * 1e9 / YEAR) for (col, rate) in st.wells]
    dh = w.subscale.active[] ? v.widgets.rugosity.value[] : 0.0
    steps(t) = round(Int, num(t))
    c = setup_injection_case(f.a, wells; fluid, boundary = copy(st.types[]),
                             temperature = T, rugosity = dh,
                             inj_time = num(w.inj) * YEAR, inj_steps = steps(w.inj_n),
                             mig_time = num(w.mig) * YEAR, mig_steps = steps(w.mig_n))
    r.run[] = (case = c, fluid = fluid, f = f, stop = Ref(false))
    r.running[] = true
    empty!(r.states)
    empty!(r.inv)
    results_surface!(r, f.s)
    r.ax.title[] = f.s === v.current[].s ? v.current[].title : ""
    r.status.text[] = "Running ..."
    background(() -> simulate!(r))
end

# The simulation of r.run[]: each report step shown as it comes (the time
# slider follows), then the trapping inventory
function simulate!(r)
    run = r.run[]
    res = try
        simulate_injection_case(run.case; info_level = -1,
            on_report = function (k, out)
                push!(r.states, out[:Reservoir])
                r.time.range[] = 1:k
                set_close_to!(r.time, k)
                show_state!(r)
                yield()
                return !run.stop[]
            end)
    finally
        r.running[] = false
    end
    if isnothing(res)
        r.status.text[] = "Stopped after $(length(r.states)) report steps."
    else
        history = trapping_history(res, run.case.model, run.fluid, run.f.s, run.f.traps)
        # the legend beside, under the colour bar
        plot_trapping_history!(r.inv, history; legend = false)
        Legend(r.side[6, 1], r.inv; halign = :left, framevisible = false)
        r.status.text[] = @sprintf("Done: %d report steps.", length(r.states))
    end
end

# --------------------------------- An app -------------------------------------

# A setup viewer (Utsira at coarsening 4, exploreSimulation's start) and its
# results figure, 850 pixels high to fit a laptop's browser window. Returns
# (v, st, r): the viewer, the setup's state, the results figure.
function new_app(; size = (1150, 850))
    name, c = "Utsirafm", 4
    f0 = formation(name, c)
    st = new_setup(f0)
    r = results_figure((900, size[2]))
    vref = Ref{Any}(nothing)
    v = interactive_traps(f0.s; res = f0.traps, charts = false, lighting = true, size,
                          title = "$name, coarsening $c",
                          widgets = g -> setup_widgets(g, name, c),
                          on_click = (col, b) -> on_click!(vref[], st, col, b),
                          on_draw = v -> overlay_plots!(v, st))
    vref[] = v
    connect_app!(v, st, r)
    return (v = v, st = st, r = r)
end

# The app's widgets' callbacks
function connect_app!(v, st, r)
    w = v.widgets
    # another formation or coarsening: its surface, no wells, all open
    for m in (w.formation, w.coarsening)
        on(m.selection) do _
            background() do
                n, k = w.formation.selection[], w.coarsening.selection[]
                v.info.text[] = "Reading $n ..."
                f = formation(n, k)
                st.f[] = f
                st.types[] = fill(:open, length(f.sides))
                empty!(st.wells)
                st.selected[] = 0
                st.pending[] = 0
                interactive_traps!(v, f.s; res = f.traps, title = "$n, coarsening $k")
            end
        end
    end
    on(w.set_all.clicks) do _
        st.types[] .= BOUNDARY_TYPES[w.btype.i_selected[]]
        update_overlays!(v, st)
    end
    on(_ -> (empty!(st.wells); st.selected[] = 0; update_overlays!(v, st)),
       w.clear.clicks)
    # the rate slider sets the selected well's rate
    on(w.rate.value) do rate
        k = st.selected[]
        k == 0 && return
        st.wells[k] = (st.wells[k][1], rate)
        update_overlays!(v, st)
    end
    # another mode: no stretch begun, the help anew
    on(_ -> (st.pending[] = 0; update_overlays!(v, st)), w.mode.selection)
    on(_ -> update_overlays!(v, st), w.btype.selection)
    on(_ -> background(() -> launch!(v, st, r)), w.run.launch.clicks)
    on(_ -> r.running[] && (r.run[].stop[] = true), w.run.stop.clicks)
    # the fringe does not combine with rugosity: one toggle switches the other off
    on(w.run.fringe.active) do on
        on && w.run.subscale.active[] && (w.run.subscale.active[] = false)
    end
    on(w.run.subscale.active) do on
        on && w.run.fringe.active[] && (w.run.fringe.active[] = false)
    end
    on(_ -> show_state!(r), r.var.selection)
    on(_ -> show_state!(r), r.time.value)
end

# ----------------------------------- Show -------------------------------------

if backend == :WGLMakie
    # an app built once before serving: compiles it, so that the first page
    # opens sooner
    new_app()

    # a local web server; every page gets an app of its own (both figures side
    # by side), live while Julia runs
    DOM = WGLMakie.Bonito.DOM
    page() = (a = new_app(); DOM.div(a.v.fig, a.r.fig; style = "display: flex"))
    server = WGLMakie.Bonito.Server(WGLMakie.Bonito.App(page), "127.0.0.1", 8080)
    println("Open http://localhost:8080 (Ctrl-C to stop)")
    wait(server)
else
    a = new_app()
    display(GLMakie.Screen(), a.r.fig)
    wait(display(GLMakie.Screen(), a.v.fig))
end

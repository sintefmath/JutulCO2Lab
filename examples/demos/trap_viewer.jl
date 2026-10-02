# A first interactive demo (roadmap Stage 11, testing the GUI approach): a
# formation of the CO2 Storage Atlas, its structural traps and spill paths, on
# the caprock in 3D. Choose the formation, the coarsening and the trap
# analysis's method; turn and zoom the caprock with the mouse; click a trap to
# see its volume, spill depth and the traps it spills into (highlighted: the
# trap red, those downstream very dark red, the spill paths between them
# thick).
#
# Written with Makie's own widgets only, so the same code runs in a window
# (GLMakie, the default) or in a browser (WGLMakie; each page gets a viewer of
# its own):
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/demos/trap_viewer.jl
#   JUTULCO2LAB_BACKEND=WGLMakie scripts/julia --project=JutulCO2Lab/examples \
#       JutulCO2Lab/examples/demos/trap_viewer.jl      (then open http://localhost:8080)
#
# Mouse: left drag turns, wheel zooms (into the scene, not cut off at the
# axis box), right drag pans; a left click (no drag) picks a trap, a right
# click clears the selection.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Printf
# in the browser, the canvas fills the window
backend == :WGLMakie && WGLMakie.activate!(resize_to = :body)
# JUTULCO2LAB_DEBUG=1 logs the events (testing the browser version)
debug(msg...) = get(ENV, "JUTULCO2LAB_DEBUG", "") == "1" &&
    (println(stderr, "[viewer] ", msg...); flush(stderr))

# the Makie extension's surface triangulation (columns of the mesh vertices)
const EXT = Base.get_extension(JutulCO2Lab, :JutulCO2LabMakieExt)

# Work started by a widget or the mouse runs as a task: in the browser a
# callback holds up the messages the plots and picking need (errors in a task
# are shown, not lost)
background(f) = @async try
    f()
catch err
    @error "trap viewer" exception = (err, catch_backtrace())
end

# Analyses done, by (formation, coarsening, method), shared by all viewers: a
# formation shown again is not read and analysed anew
const ANALYSES = Dict{Tuple{String, Int, Symbol}, Any}()
const ANALYSES_LOCK = ReentrantLock()

# ---------------------------------- A viewer ----------------------------------

# A new viewer: its figure, the 3D axis, the widgets, and its state (the
# current formation's analysis, the highlight's plots, where a button went
# down). Returns a NamedTuple.
function new_viewer()
    # no padding: the 3D scene reaches the window's edges
    fig = Figure(size = (1700, 1000), figure_padding = 0)
    # :free zooms the whole scene, so close-ups are not clipped at the axis box
    ax = Axis3(fig[1, 1]; zreversed = true, viewmode = :free, clip = false,
               xlabel = "x (m)", ylabel = "y (m)", zlabel = "depth (m)",
               xticks = LinearTicks(3), yticks = LinearTicks(3))
    # the panel floats over the axis's right side: a 3D scene is cut off at the
    # edge of its axis, so the axis takes the whole window
    panel = GridLayout(fig[1, 1]; width = 300, halign = :right, valign = :top,
                       tellwidth = false, tellheight = false,
                       alignmode = Outside(15))
    # a backdrop under the widgets (made first, so drawn below them)
    Box(panel[1:10, 1]; color = (:white, 0.85), strokevisible = false,
        alignmode = Outside(-10))

    # plain dropdowns (Makie's menus are by default also search boxes)
    menu(pos; kw...) = Menu(pos; searchable = false, kw...)
    Label(panel[1, 1], "Formation"; halign = :left)
    formation = menu(panel[2, 1]; options = atlas_grid_formations(),
                     default = "Johansenfm")
    Label(panel[3, 1], "Coarsening"; halign = :left)
    coarsening = menu(panel[4, 1]; options = 1:6, default = 2)  # (the 2nd: 2)
    Label(panel[5, 1], "Trap analysis"; halign = :left)
    method = menu(panel[6, 1]; options = ["node", "cell"], default = "node")
    Label(panel[7, 1], "Vertical exaggeration"; halign = :left)
    # (an integer default is a position in the options: here 10)
    exaggeration = menu(panel[8, 1]; options = [1, 2, 5, 10, 20, 50], default = 4)
    summary = Label(panel[9, 1], ""; halign = :left, justification = :left)
    info = Label(panel[10, 1], "Click a trap."; halign = :left, justification = :left)
    # the title floats too (an axis title would take room from the scene)
    title = Label(fig[1, 1], ""; font = :bold, valign = :top, padding = (0, 0, 0, 15),
                  tellwidth = false, tellheight = false)
    # in front of the 3D scene (2D blocks are otherwise hidden where it is
    # nearer in the depth buffer)
    for b in [title; contents(panel)]
        translate!(b.blockscene, 0, 0, 1000)
    end

    v = (fig = fig, ax = ax, title = title, formation = formation, coarsening = coarsening,
         method = method, exaggeration = exaggeration, summary = summary,
         info = info, current = Ref{Any}(nothing), highlight = Any[],
         pressed_at = Dict{Mouse.Button, Point2f}())
    connect!(v)
    analyse!(v)
    return v
end

# The widgets' and the mouse's callbacks
function connect!(v)
    for menu in (v.formation, v.coarsening, v.method)
        on(_ -> background(() -> analyse!(v)), menu.selection)
    end
    on(v.exaggeration.selection) do _
        isnothing(v.current[]) || set_aspect!(v, v.current[].s)
    end
    on(ev -> mouse!(v, ev), events(v.ax.scene).mousebutton)
end

# ------------------------------ The trap analysis -----------------------------

# The axis to scale in x and y, its depth exaggerated
function set_aspect!(v, s)
    span(x) = -(reverse(extrema(x))...)
    m = max(span(s.x), span(s.y))
    v.ax.aspect[] = (span(s.x) / m, span(s.y) / m,
                     v.exaggeration.selection[] * max(span(s.z), 1.0) / m)
end

# The trap surface, trap analysis and trap volumes of a formation, from the
# cache if done before: a Tuple (shared, not to be modified)
function analysis(name, c, m)
    lock(ANALYSES_LOCK) do
        get!(ANALYSES, (name, c, m)) do
            mesh = atlas_mesh(name; coarsening = c)
            s = trap_surface(mesh.mesh, mesh.columns)
            res = trap_analysis(s; method = m)
            (s, res, trap_volumes(s, res))
        end
    end
end

# The selected formation's traps, drawn; the analysis and the caprock mesh plot
# with its vertices' columns (for picking) kept in `v.current`
function analyse!(v)
    name, c = v.formation.selection[], v.coarsening.selection[]
    m = Symbol(v.method.selection[])
    debug("analysing ", name, ", coarsening ", c, ", ", m)
    v.summary.text[] = "Computing $name ..."
    s, res, vol = analysis(name, c, m)

    # redraw: spill regions, traps, spill paths
    empty!(v.ax)
    empty!(v.highlight)
    plot_traps!(v.ax, s, res)
    surface = first(p for p in v.ax.scene.plots if p isa Mesh)
    _, _, vcols = EXT._surface_triangles(s, 0.0)
    v.current[] = (s = s, res = res, vol = vol, surface = surface, vcols = vcols)

    set_aspect!(v, s)
    v.title.text[] = "$name, coarsening $c, $(m)-based traps"
    reset_limits!(v.ax)
    v.summary.text[] = @sprintf("%d columns\n%d traps, %.3g Mm³ in all",
                                length(s.zc), length(vol), sum(vol) / 1e6)
    v.info.text[] = "Click a trap."
    debug("drawn: ", length(s.zc), " columns")
end

# ------------------------------ Selecting a trap ------------------------------

function describe(c, t)
    down = downstream_traps(c.res, t)[2:end]
    txt = @sprintf("Trap %d\nvolume %.3g Mm³\nspill depth %.0f m", t, c.vol[t] / 1e6,
                   c.res.trap_z[t])
    isempty(down) && return txt * "\nspills out of the domain"
    return txt * @sprintf("\nspills into %d trap(s)\nholding %.3g Mm³ more",
                          length(down), sum(c.vol[down]) / 1e6)
end

# Trap `t` filled red, the traps downstream very dark red, the spill paths
# from them thick, all lifted slightly above the caprock (t = 0: none)
function highlight!(v, t)
    foreach(p -> delete!(v.ax, p), v.highlight)
    empty!(v.highlight)
    t == 0 && return
    s, res = v.current[].s, v.current[].res
    δ = EXT._lift(s)
    traps = downstream_traps(res, t)
    for (ts, colour) in ((traps[1:1], :red), (traps[2:end], RGBf(0.35, 0, 0)))
        cols = findall(in(Set(ts)), res.traps)
        isempty(cols) && continue
        sub = (; s..., polygons = s.polygons[cols])
        plot_ve_surface!(v.ax, sub; color = colour, offset = -δ, strokewidth = 0)
        push!(v.highlight, v.ax.scene.plots[end])
    end
    pts = [Point3f(s.cx[k], s.cy[k], s.zc[k] - 2δ) for u in traps for r in res.rivers[u]
           for (a, b) in zip(r[1:end-1], r[2:end]) for k in (a, b)]
    isempty(pts) ||
        push!(v.highlight, linesegments!(v.ax, pts; color = :red3, linewidth = 4))
end

# A click is a press and release without moving (a drag turns or pans the
# view): left picks a trap, right clears the selection
function mouse!(v, ev)
    ev.button in (Mouse.left, Mouse.right) || return
    pos = Point2f(events(v.ax.scene).mouseposition[])
    debug(ev.button, " ", ev.action, " at ", pos)
    if ev.action == Mouse.press
        v.pressed_at[ev.button] = pos
        return
    end
    sum(abs2, pos .- get(v.pressed_at, ev.button, Point2f(Inf))) < 9 || return
    isnothing(v.current[]) && return
    if ev.button == Mouse.right
        highlight!(v, 0)
        v.info.text[] = "Click a trap."
        return
    end
    c = v.current[]
    background(() -> select!(v, c, pos))
end

# The trap under screen position `pos`: picked on the caprock (in the browser
# picking asks the page, so this runs as a task)
function select!(v, c, pos)
    plot, k = pick(v.ax.scene, pos)
    debug("picked at ", pos, ": ", plot === c.surface ? "caprock" : string(plot), k)
    # only the caprock itself identifies a column
    plot === c.surface || return
    col = c.vcols[k]
    t = c.res.traps[col]
    if t == 0
        r = c.res.trap_regions[col]
        v.info.text[] = r == 0 ? "No trap here; spills out of the domain." :
                                 "No trap here; spills into trap $r."
        highlight!(v, r)
        return
    end
    v.info.text[] = describe(c, t)
    highlight!(v, t)
end

# ----------------------------------- Show -------------------------------------

if backend == :WGLMakie
    # a viewer built once before serving: compiles the plotting, caches the
    # first formation, so that the first page opens sooner
    new_viewer()

    # a local web server; every page gets a viewer of its own, live while
    # Julia runs
    app = WGLMakie.Bonito.App(() -> new_viewer().fig)
    server = WGLMakie.Bonito.Server(app, "127.0.0.1", 8080)
    println("Open http://localhost:8080 (Ctrl-C to stop)")
    wait(server)
else
    screen = display(new_viewer().fig)
    wait(screen)
end

# A first interactive demo (roadmap Stage 11): the formations of the CO2
# Storage Atlas in the interactive trap viewer, `interactive_traps` (co2lab's
# interactiveTrapping), with menus for the formation and its coarsening.
# Turn and zoom the caprock with the mouse; click a column to see the traps
# CO2 injected there fills (shift-click: the traps spilling into it), their
# volumes in a pie and a bar chart; the panel's switches show or hide traps,
# spill regions (off: depth colours), lighting, contours and the red path.
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
# axis box), right drag pans; a left click (no drag) picks a column, a right
# click clears the selection.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
# in the browser, the canvas fills the window
backend == :WGLMakie && WGLMakie.activate!(resize_to = :body)

# Grids read, by (formation, coarsening), shared by all viewers: a formation
# shown again is not read anew
const SURFACES = Dict{Tuple{String, Int}, Any}()
const SURFACES_LOCK = ReentrantLock()

# The trap surface of an Atlas formation, from the cache if read before
function surface(name, c)
    lock(SURFACES_LOCK) do
        get!(SURFACES, (name, c)) do
            mesh = atlas_mesh(name; coarsening = c)
            trap_surface(mesh.mesh, mesh.columns)
        end
    end
end

# Reading a grid runs as a task: in the browser a callback holds up the
# messages the plots need (errors in a task are shown, not lost)
background(f) = @async try
    f()
catch err
    @error "trap viewer" exception = (err, catch_backtrace())
end

# ---------------------------------- A viewer ----------------------------------

# The formation and coarsening menus, in layout g (the top of the viewer's
# panel). Returns the two menus.
function atlas_menus(g, name, c)
    menu(pos, options, default, w) = Menu(pos; options, default, searchable = false,
                                          width = w)
    Label(g[1, 1], "Formation"; halign = :left)
    formation = menu(g[2, 1], atlas_grid_formations(), name, 190)
    Label(g[1, 2], "Coarsening"; halign = :left)
    # (an integer default is a position in the options: here 2)
    coarsening = menu(g[2, 2], 1:6, c, 100)
    return (formation, coarsening)
end

# A viewer of Johansen at coarsening 2, with the formation and coarsening
# menus at the top of its panel. Returns the viewer (interactive_traps's
# NamedTuple).
function new_viewer()
    name, c = "Johansenfm", 2
    v = interactive_traps(surface(name, c); spill_regions = true, lighting = true,
                          title = "$name, coarsening $c",
                          widgets = g -> atlas_menus(g, name, c))
    formation, coarsening = v.widgets

    # another formation or coarsening: its surface in the same viewer
    for m in (formation, coarsening)
        on(m.selection) do _
            background() do
                n, k = formation.selection[], coarsening.selection[]
                v.info.text[] = "Reading $n ..."
                interactive_traps!(v, surface(n, k); title = "$n, coarsening $k")
            end
        end
    end
    return v
end

# ----------------------------------- Show -------------------------------------

if backend == :WGLMakie
    # a viewer built once before serving and put through its paces (a column
    # selected, each switch, the trap detail, the other method): compiles
    # them, so that the first page and its first uses are quick
    v = new_viewer()
    c = v.current[]
    interactive_traps!(v, c.s; title = c.title, injection_column = c.res.top[1])
    for t in v.toggles, _ in 1:2
        t.active[] = !t.active[]
    end
    foreach(_ -> v.detail.clicks[] += 1, 1:2)
    for k in (2, 1)
        v.method.i_selected[] = k
        sleep(0.5)
        while startswith(v.info.text[], "Computing")
            sleep(0.5)
        end
    end

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

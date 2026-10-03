# co2lab-demo's exploreCapacity (roadmap Stage 11 item 2): the CO2 storage
# capacity of a formation of the CO2 Storage Atlas, column by column, as it
# changes with the conditions — the temperature gradient, the pressure's
# deviation from hydrostatic, the residual saturations and the dissolution —
# shown on the caprock in the interactive trap viewer (`interactive_traps`).
#
# The estimate (`trapping_capacity`) is an upper bound: every structural trap
# filled to its spill depth, all the brine saturated with dissolved CO2. In
# the traps: free CO2 and CO2 dissolved in the residual brine; below them:
# residual CO2 and CO2 dissolved in the brine. The CO2 density is taken at
# the caprock's pressure and temperature (JutulDarcy's tables).
#
# Choose the formation, its coarsening and what colours the caprock: the
# trap analysis, or a map — caprock depth, thickness, temperature, pressure,
# CO2 density, total capacity per area, each trap's structural capacity, the
# structural capacity reachable from each column (its trap and those
# downstream).
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/demos/explore_capacity.jl
#   JUTULCO2LAB_BACKEND=WGLMakie scripts/julia --project=JutulCO2Lab/examples \
#       JutulCO2Lab/examples/demos/explore_capacity.jl   (then open http://localhost:8080)

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Printf
# in the browser, the canvas fills the window
backend == :WGLMakie && WGLMakie.activate!(resize_to = :body)

# Reading a grid runs as a task: in the browser a callback holds up the
# messages the plots need (errors in a task are shown, not lost)
background(f) = @async try
    f()
catch err
    @error "explore capacity" exception = (err, catch_backtrace())
end

# ------------------------------- The formations -------------------------------

# Grids read, by (formation, coarsening), shared by all viewers: the trap
# surface and the porosity per column
const FORMATIONS = Dict{Tuple{String, Int}, Any}()
const FORMATIONS_LOCK = ReentrantLock()

# A formation's trap surface and porosity, from the cache if read before; a
# formation without Atlas rock values gets porosity 0.2 (co2lab's default)
function formation(name, c)
    lock(FORMATIONS_LOCK) do
        get!(FORMATIONS, (name, c)) do
            a = try
                atlas_ve_domain(name; coarsening = c)
            catch err
                err isa ArgumentError || rethrow()
                atlas_ve_domain(name; coarsening = c, poro = 0.2, perm = 1e-13)
            end
            (s = trap_surface(a.mesh, a.columns), poro = a.ve.domain[:porosity])
        end
    end
end

# ------------------------------ The estimate ----------------------------------

# exploreCapacity's conditions not on a slider: the sea floor and the water
const SEAFLOOR_DEPTH, SEAFLOOR_TEMP, WATER_DENSITY = 100.0, 7.0, 1000.0

# The sliders: label, range, start value (exploreCapacity's)
const SLIDERS = [
    (label = "T gradient (°C/km)", range = 10:0.2:50, startvalue = 35.6),
    # deviation from hydrostatic pressure
    (label = "Pressure dev. (%)", range = -50:1:50, startvalue = 0),
    (label = "Residual CO2", range = 0:0.01:0.5, startvalue = 0.21),
    (label = "Residual brine", range = 0:0.01:0.5, startvalue = 0.11),
    # dissolved CO2 per brine volume when saturated
    (label = "Dissolution (kg/m³)", range = 0:1:100, startvalue = 53)]

# The caprock's pressure (Pa; hydrostatic from the sea surface, no
# atmosphere, as co2lab), temperature (K) and the capacity, for the sliders'
# values `x`
function estimate(s, res, poro, x)
    grad, dev, sg, sw, dis = x
    p = s.zc .* WATER_DENSITY .* 9.80665 .* (1 + dev / 100)
    T = SEAFLOOR_TEMP .+ (s.zc .- SEAFLOOR_DEPTH) .* grad ./ 1000 .+ 273.15
    cap = trapping_capacity(s, res; p, T, res_gas = sg, res_water = sw,
                            dis_max = dis, poro)
    return (p = p, T = T, cap = cap)
end

# The maps: their names and values per column (nothing: the trap analysis)
const MAPS = ["caprock topography (traps)", "caprock depth (m)",
              "formation thickness (m)", "caprock temperature (°C)",
              "caprock pressure (MPa)", "CO2 density (kg/m³)",
              "total capacity (t/m²)", "structural trap capacity (Mt)",
              "reachable structural capacity (Mt)"]

function map_values(k, s, res, e)
    k == 1 && return nothing
    k == 2 && return s.zc
    k == 3 && return s.H
    k == 4 && return e.T .- 273.15
    k == 5 && return e.p ./ 1e6
    k == 6 && return e.cap.density
    total = e.cap.structural .+ e.cap.residual .+ e.cap.dissolved
    k == 7 && return total ./ s.area ./ 1e3
    k == 9 && return reachable_capacity(res, e.cap.structural) ./ 1e9
    # each trap's structural capacity on its columns, NaN elsewhere
    per_trap = zeros(length(res.trap_z))
    for (c, t) in enumerate(res.traps)
        t > 0 && (per_trap[t] += e.cap.structural[c])
    end
    return [t == 0 ? NaN : per_trap[t] / 1e9 for t in res.traps]
end

# exploreCapacity's report: the total and its parts
function report(cap)
    parts = (sum(cap.structural), sum(cap.residual), sum(cap.dissolved))
    total = sum(parts)
    line(name, x) = @sprintf("%s: %.2f Gt (%.1f%%)", name, x / 1e12, 100x / total)
    return @sprintf("Total trapping capacity: %.2f Gt\n", total / 1e12) *
           join(line.(("Structural", "Residual", "Dissolved"), parts), "\n")
end

# ---------------------------------- A viewer ----------------------------------

# The app's own widgets, in layout g (the top of the viewer's panel): the
# formation and coarsening, the map, the sliders, the report. Returns them.
function capacity_widgets(g, name, c)
    menu(pos, options, default, w) = Menu(pos; options, default, searchable = false,
                                          width = w)
    Label(g[1, 1], "Formation"; halign = :left)
    f = menu(g[2, 1], atlas_grid_formations(), name, 190)
    Label(g[1, 2], "Coarsening"; halign = :left)
    # (an integer default is a position in the options: here c)
    k = menu(g[2, 2], 1:6, c, 100)
    Label(g[3, 1:2], "Colour by"; halign = :left)
    m = menu(g[4, 1:2], MAPS, MAPS[8], 300)
    sliders = SliderGrid(g[5, 1:2], SLIDERS...; tellwidth = false)
    text = Label(g[6, 1:2], ""; halign = :left, justification = :left)
    return (formation = f, coarsening = k, map = m, sliders = sliders, report = text)
end

# Utsira at coarsening 2 (exploreCapacity's start), its structural traps'
# capacity shown. Returns the viewer (interactive_traps's NamedTuple).
function new_viewer()
    name, c = "Utsirafm", 2
    v = interactive_traps(formation(name, c).s; spill_regions = true, lighting = true,
                          charts = false, selection = false,
                          title = "$name, coarsening $c",
                          widgets = g -> capacity_widgets(g, name, c))
    w = v.widgets
    # the formation on screen (its porosity goes with its surface)
    shown = Ref(formation(name, c))
    refresh() = refresh!(v, shown[])
    # another formation or coarsening: its surface in the same viewer
    for m in (w.formation, w.coarsening)
        on(m.selection) do _
            background() do
                n, k = w.formation.selection[], w.coarsening.selection[]
                v.info.text[] = "Reading $n ..."
                f = formation(n, k)
                interactive_traps!(v, f.s; title = "$n, coarsening $k")
                shown[] = f
                refresh()
            end
        end
    end
    # another map, other conditions: the estimate anew
    on(_ -> refresh(), w.map.selection)
    foreach(sl -> on(_ -> refresh(), sl.value), w.sliders.sliders)
    refresh()
    return v
end

# The estimate for the formation `f` on screen and the sliders: the report,
# the map
function refresh!(v, f)
    w = v.widgets
    f.s === v.current[].s || return
    s, res = v.current[].s, v.current[].res
    e = estimate(s, res, f.poro, [sl.value[] for sl in w.sliders.sliders])
    w.report.text[] = report(e.cap)
    k = findfirst(==(w.map.selection[]), MAPS)
    show_column_values!(v, map_values(k, s, res, e); label = MAPS[k])
end

# ----------------------------------- Show -------------------------------------

if backend == :WGLMakie
    # a viewer built once before serving: compiles it, so that the first page
    # opens sooner
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

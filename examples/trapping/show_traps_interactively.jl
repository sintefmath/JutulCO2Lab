# co2lab-spillpoint's showTrapsInteractively: the interactive trap viewer
# (`interactive_traps`, co2lab's interactiveTrapping) on three caprocks — a
# small synthetic one, the Johansen formation of the CO2 Storage Atlas and an
# IGEMS surface.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/trapping/show_traps_interactively.jl
#
# In the viewer (GLMakie, the default; close a window to go on to the next):
# - left click on the caprock: the trap CO2 injected there rises to, and the
#   traps it spills on into, coloured along the path; a pie chart of the
#   volume in the first trap, along the path and in the other traps, and the
#   path traps' volumes as bars on a log scale — far-off traps are reached
#   by a slower front and rest on more uncertain data;
# - shift-click or middle click: the traps spilling into it instead (possible
#   injection sites that feed a large trap); right click: clear;
# - the panel: the method, the vertical exaggeration, traps, spill regions
#   (off: depth colours), lighting, depth contours, the red path, a view from
#   above (the picking works best from above), the starting view, and the
#   selected trap alone (its fill red, its ceiling blue).
# JUTULCO2LAB_BACKEND=CairoMakie writes the viewers' first views to files.
# co2lab's "Setup VE simulation" button is roadmap Stage 11 item 3.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend

# co2lab's choice (useCell = true); its viewer runs the node-based method
# whatever it is asked (co2lab-issues.org), ours the one asked for
method = :cell

# ------------------------------ 1. A synthetic grid ---------------------------

# 101 x 101 cells over 1 km x 1 km, 10 m thick, at 100 m depth: dipping in x
# and perturbed by a periodic function of x and y, which makes a hierarchy of
# small traps (co2lab: a cartGrid with its nodes moved)
L, Lp, H = 1000.0, 200.0, 10.0
x = y = range(0, L; length = 102)
top = [100 + 0.01a - 2sin(pi * a / Lp) * sin(pi * b / Lp) for a in x, b in y]
grid = surface_mesh(x, y, top; bottom = top .+ H)
s = trap_surface(grid.mesh, grid.columns)
v = interactive_traps(s; method, lighting = true, view = (180, 50),
                      title = "A synthetic caprock")
show_figure(v.fig, "interactive_synthetic")

# --------------------------- 2. The Johansen formation ------------------------

# The Atlas grid, coarsened twice (the full resolution is more accurate), with
# depth contours
johansen = atlas_mesh("Johansenfm"; coarsening = 2)
s = trap_surface(johansen.mesh, johansen.columns)
v = interactive_traps(s; method, lighting = true, contours = true, view = (-70, 80),
                      title = "Johansen formation (CO2 Storage Atlas, coarsening 2)")
show_figure(v.fig, "interactive_johansen")

# ------------------------------ 3. An IGEMS surface ---------------------------

# IGEMS: alternative caprock morphologies, made stochastically from
# stratigraphic and structural scenarios (http://files.nr.no/igems/data.pdf);
# surface 1 of the OSS family, every 2nd value in each direction
igems = igems_mesh("OSS", 1; coarsening = 2)
s = trap_surface(igems.mesh, igems.columns)
v = interactive_traps(s; method, view = (180, 50), title = "IGEMS surface OSS 1")
show_figure(v.fig, "interactive_igems")

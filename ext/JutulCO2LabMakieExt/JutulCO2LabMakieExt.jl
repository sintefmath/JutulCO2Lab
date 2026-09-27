"""
Makie extension of JutulCO2Lab: map views of VE domains and CO2 Storage Atlas
data (roadmap Stage 5 item 5), and of trapping analyses (Stage 8). Loaded
with any Makie backend (GLMakie, CairoMakie). The functions are documented
where they are declared, in `src/plotting.jl`.
"""
module JutulCO2LabMakieExt

using JutulCO2Lab, Makie
import JutulCO2Lab: plot_ve_map, plot_ve_map!, plot_atlas_map, plot_atlas_map!,
    plot_atlas_context!, plot_ve_outline!, plot_ve_surface!, show_figure,
    plot_traps, plot_traps!, ve_axis3, plot_trapping_history, plot_trapping_history!

include("map_plots.jl")
include("atlas_plots.jl")
include("trap_plots.jl")
include("surface_3d.jl")
include("inventory_plots.jl")

# ------------------------------ Figures, backends -----------------------------

# A figure with one map axis (equal scaling of x and y)
function _map_figure(title; size = (900, 800))
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title = title, aspect = DataAspect(), xlabel = "x (m)",
              ylabel = "y (m)", xticks = LinearTicks(4))
    return (fig, ax)
end

function show_figure(fig, name; dir = "figures")
    if nameof(Makie.current_backend()) == :CairoMakie
        mkpath(dir)
        path = joinpath(dir, "$name.png")
        save(path, fig)
        println("saved $path")
    elseif nameof(Makie.current_backend()) == :GLMakie
        # a window of its own; a script's windows close when Julia exits, so
        # outside the REPL wait until this one is closed
        screen = Makie.current_backend().Screen()
        display(screen, fig)
        isinteractive() || wait(screen)
    else
        display(fig)
    end
    return nothing
end

end # module

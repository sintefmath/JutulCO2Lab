"""
Makie extension of JutulCO2Lab: map views of VE domains and CO2 Storage Atlas
data (roadmap Stage 5 item 5). Loaded with any Makie backend (GLMakie,
CairoMakie). The functions are documented where they are declared, in
`src/plotting.jl`.
"""
module JutulCO2LabMakieExt

using JutulCO2Lab, Makie
import JutulCO2Lab: plot_ve_map, plot_ve_map!, plot_atlas_map, plot_atlas_map!,
    plot_atlas_context!, plot_ve_outline!, show_figure

include("map_plots.jl")
include("atlas_plots.jl")

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
    else
        display(fig)
    end
    return nothing
end

end # module

# Map view of a trapping analysis (co2lab's trap plots in its spill-point
# examples, without the interactive viewer).

function plot_traps(s, res; title = "", colorbar = true, kwarg...)
    fig, ax = _map_figure(title)
    plot_traps!(ax, s, res; kwarg...)
    p = _background_plot(ax)
    colorbar && !isnothing(p) && Colorbar(fig[1, 2], p)
    return fig
end

function plot_traps!(ax, s, res; regions = true, trap_values = nothing,
                     colormap = :viridis, rivers = true, labels = false,
                     outline = true, faults = true)
    polys = _polygons(s)
    palette = Makie.to_colormap(:tab20)
    colour(k) = palette[mod1(k, length(palette))]
    light(c) = RGBAf(0.45c.r + 0.55, 0.45c.g + 0.55, 0.45c.b + 0.55, 1)

    # spill regions in light colours, the region draining out in grey
    if regions
        grey = RGBAf(0.9, 0.9, 0.9, 1)
        rc = [k == 0 ? grey : light(colour(k)) for k in res.trap_regions]
        poly!(ax, polys; color = rc, strokewidth = 0)
    end

    # traps: coloured by `trap_values` (one per trap), else as their regions
    inside = findall(>(0), res.traps)
    if isnothing(trap_values)
        poly!(ax, polys[inside]; color = colour.(res.traps[inside]), strokewidth = 0)
    else
        values = Float64[trap_values[res.traps[c]] for c in inside]
        poly!(ax, polys[inside]; color = values, colormap = colormap, strokewidth = 0)
    end

    # rivers through the column centroids
    if rivers
        pts = [Point2f(s.cx[c], s.cy[c]) for rt in res.rivers for r in rt
               for (a, b) in zip(r[1:end-1], r[2:end]) for c in (a, b)]
        linesegments!(ax, pts; color = :black, linewidth = 1)
    end
    if labels
        text!(ax, s.cx[res.top], s.cy[res.top];
              text = string.(res.traps[res.top]), align = (:center, :center),
              fontsize = 11)
    end
    outline && plot_ve_outline!(ax, s)
    # sealing faults over the outline
    if faults && hasproperty(s, :sealed) && any(s.sealed)
        pts = [Point2f(s.x[n], s.y[n]) for k in findall(s.sealed) for n in s.edges[:, k]]
        linesegments!(ax, pts; color = :red3, linewidth = 3)
    end
    return ax
end

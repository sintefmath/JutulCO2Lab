# Vertical cross-sections of VE and hybrid grids (ve_section): fills by block
# type, segments by kind.

# The segment kinds, drawn in this order (barriers on top), with their legend
# labels (nothing: none)
const _SECTION_KINDS = ((:internal, nothing), (:outline, nothing),
                        (:leaky, "leaky barrier"), (:sealed, "sealed barrier"))

# plot_ve_section: documented in src/plotting.jl
function plot_ve_section(section; exaggeration = 10, width = 1200, title = "",
                         legend = true, kwarg...)
    s = [p[1] for poly in section.polygons for p in poly]
    z = [p[2] for poly in section.polygons for p in poly]
    # the initial height that gives the exaggeration at this width (the margins
    # roughly allowed for), within reason
    span(x) = max(maximum(x) - minimum(x), 1.0)
    span_s, span_z = span(s), span(z)
    height = clamp(round(Int, (width - 150) * exaggeration * span_z / span_s) + 150,
                   250, 1000)
    fig = Figure(size = (width, height))
    ax = Axis(fig[1, 1]; title, xlabel = "distance along the section (m)",
              ylabel = "depth (m)", yreversed = true)
    plot_ve_section!(ax, section; kwarg...)
    # the legend beside the axis, not over the section
    legend && Legend(fig[1, 2], ax)
    return fig
end

# plot_ve_section!: documented in src/plotting.jl
function plot_ve_section!(ax, section; ve_color = (:steelblue, 0.35),
                          fine_color = (:darkorange, 0.5),
                          outline = (color = :black, linewidth = 1),
                          internal = (color = :grey60, linewidth = 0.5),
                          sealed = (color = :red, linewidth = 2.5),
                          leaky = (color = :red, linewidth = 2.5, linestyle = :dash))
    # the fills: VE blocks' cells, fine 3D cells
    polys = [Point2d.(collect(p)) for p in section.polygons]
    for (ve, color, label) in ((true, ve_color, "VE"), (false, fine_color, "3D"))
        sel = findall(==(ve), section.is_ve)
        isempty(sel) || poly!(ax, polys[sel]; color, strokewidth = 0, label)
    end

    # the segments by kind; a kind whose style is nothing is left out
    styles = (; outline, internal, sealed, leaky)
    for (kind, label) in _SECTION_KINDS
        style = styles[kind]
        sel = findall(==(kind), section.kind)
        (isnothing(style) || isempty(sel)) && continue
        segs = [Point2d(a) => Point2d(b) for (a, b) in section.segments[sel]]
        linesegments!(ax, segs; label, style...)
    end
    return ax
end

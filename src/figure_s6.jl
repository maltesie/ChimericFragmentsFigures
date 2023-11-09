function plot_figure_s6(assets_folder::String, interact::InteractionsNew)
    resfactor = 1.

    fig_1 = Figure(resolution=(900*resfactor, 1200*resfactor))

    ax_graph = Axis(fig_1[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(assets_folder, "IGR_CDS_graph_rot90.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())

    Label(fig_1[1,1, TopLeft()], "a", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save( "figure_s6_1.svg", fig_1)
    save( "figure_s6_1.png", fig_1, px_per_unit = 3)

    fig_2 = Figure(resolution=(1200*resfactor, 1200*resfactor))

    colors = Makie.wong_colors()

    ax_graph = Axis(fig_2[1,1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(assets_folder, "hfq_lcd_igr_vcr069_vc1803.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())


    ga = fig_2[2:3,1] = GridLayout(alignmode=Outside())
    ax1_1 = Axis(ga[1, 1], ylabel="count", title="basepairing predictions")

    limits_low, limits_high = 1530, 1630

    aggregation_plot!(ax1_1, interact, "Vcr069:VC1803", 0.9)
    axislegend(ax1_1)
    hidexdecorations!(ax1_1, grid = false)
    hidespines!(ax1_1, :l, :t, :r)
    xlims!(ax1_1, limits_low, limits_high)

    idx = findfirst(interact.nodes.name .== "Vcr069:VC1803")
    ref = interact.nodes.ref[idx]
    st = interact.nodes.strand[idx]
    interval = Interval(ref, interact.nodes.left[idx], interact.nodes.right[idx], st)
    coverage_tex = Coverage(joinpath(assets_folder, "..", "figure_6", "trimmed_tex_01_1_forward.bw"), joinpath(assets_folder, "..", "figure_6", "trimmed_tex_01_1_reverse.bw"))
    coverage_term = Coverage(joinpath(assets_folder, "..", "figure_6", "trimmed_12898-mock-rep1_S33_R1_001_forward.bw"), joinpath(assets_folder, "..", "figure_6", "trimmed_12898-mock-rep1_S33_R1_001_reverse.bw"))

    ax1_2 = Axis(ga[2,1], ylabel="count", title="dRNA-seq")
    coverage_plot!(ax1_2, coverage_tex, interval, colors[1])
    hidexdecorations!(ax1_2, grid = false)
    hidespines!(ax1_2, :l, :t, :r)
    xlims!(ax1_2, limits_low, limits_high)

    ax1_3 = Axis(ga[3,1], ylabel="count", xlabel="position in IGR", title="TERM-seq")
    coverage_plot!(ax1_3, coverage_term, interval, colors[1])
    hidespines!(ax1_3, :l, :t, :r)
    xlims!(ax1_3, limits_low, limits_high)

    linkxaxes!(ax1_1, ax1_2, ax1_3)

    img4 = rotr90(load(joinpath(assets_folder, "northern_2.png")))
    ax4 = Axis(fig_2[1:3, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    hidedecorations!(ax4)
    image!(ax4, img4, aspect = DataAspect())

    Label(fig_2[1,1, TopLeft()], "b", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig_2[2,1, TopLeft()], "c", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig_2[1,2, TopLeft()], "d", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_s6_2.svg", fig_2)
    save("figure_s6_2.png", fig_2, px_per_unit = 2)
end
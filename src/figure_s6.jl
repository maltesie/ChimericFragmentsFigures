function plot_figure_s6(assets_folder::String)
    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 1600*resfactor))

    ax_graph = Axis(fig[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(assets_folder, "screenshot_graph_cropped.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())

    Label(fig[1,1, TopLeft()], "a", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    ax_table = Axis(fig[2, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_tabl = rotr90(load(joinpath(assets_folder, "screenshot_table_cropped.png")))
    hidedecorations!(ax_table)
    image!(ax_table, img_tabl, aspect = DataAspect())

    Label(fig[2,1, TopLeft()], "b", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_s6_1.svg", fig)
    save("figure_s6_1.png", fig)

    cp(joinpath(assets_folder, "figure_si_S6_2.svg"), "figure_s6_2.svg"; force=true)
    cp(joinpath(assets_folder, "figure_si_S6_2.png"), "figure_s6_2.png"; force=true)

    fig5 = Figure(resolution=(1200*resfactor, 800*resfactor))

    ax_summary = Axis(fig5[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_summary = rotr90(load(joinpath(assets_folder, "screenshot_summary_cropped_annot.png")))
    hidedecorations!(ax_summary)
    image!(ax_summary, img_summary, aspect = DataAspect())

    Label(fig5[1,1, TopLeft()], "g", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_s6_3.svg", fig5)
    save("figure_s6_3.png", fig5, px_per_unit = 2)
end

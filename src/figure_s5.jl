function plot_figure_s5(assets_folder::String)
    resfactor = 1.
    fig = Figure(resolution=(1400*resfactor, 900*resfactor))

    ax_graph = Axis(fig[1:4, 3], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(assets_folder, "screenshot_aggregation.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())

    ax_graph2 = Axis(fig[1:4, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph2 = rotr90(load(joinpath(assets_folder, "screenshot_ligation_points.png")))
    hidedecorations!(ax_graph2)
    image!(ax_graph2, img_graph2, aspect = DataAspect())

    ax_controls = Axis(fig[1:4, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_controls = rotr90(load(joinpath(assets_folder, "controls_annot.png")))
    hidedecorations!(ax_controls)
    image!(ax_controls, img_controls, aspect = DataAspect())

    ax_menu = Axis(fig[5, 1:2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_menu = rotr90(load(joinpath(assets_folder, "screenshot_plots_menu_annot.png")))
    hidedecorations!(ax_menu)
    image!(ax_menu, img_menu, aspect = DataAspect())


    Label(fig[1,1, TopLeft()], "A", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[1,2, TopLeft()], "B", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[1,3, TopLeft()], "C", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[5,1, TopLeft()], "D", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_s5.svg", fig)
    save("figure_s5.png", fig, px_per_unit = 2)
end
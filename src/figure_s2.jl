using CairoMakie, FileIO

function plot_figure_s2(assets_folder::String)
    resfactor = 1.0
    fig = Figure(resolution=(900*resfactor, 1200*resfactor))

    ax_graph = Axis(fig[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(assets_folder, "IGR_CDS_graph_rot90.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())

    Label(fig[1,1, TopLeft()], "A", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save( "figure_s2.svg", fig)
end

using CairoMakie, FileIO

function si_plot()
    resfactor = 1.0
    fig = Figure(resolution=(900*resfactor, 1200*resfactor))

    ax_graph = Axis(fig[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(@__DIR__, "IGR_CDS_graph_rot90.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())

    Label(fig[1,1, TopLeft()], "A", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save(joinpath(@__DIR__, "figure_si_S2.svg"), fig)

end
si_plot()

using CairoMakie, FileIO

function si_plot()
    resfactor = 1.0
    fig = Figure(resolution=(1200*resfactor, 1600*resfactor))

    ax_graph = Axis(fig[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(@__DIR__, "screenshot_graph_cropped.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())

    Label(fig[1,1, TopLeft()], "A", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    ax_table = Axis(fig[2, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_tabl = rotr90(load(joinpath(@__DIR__, "screenshot_table_cropped.png")))
    hidedecorations!(ax_table)
    image!(ax_table, img_tabl, aspect = DataAspect())

    Label(fig[2,1, TopLeft()], "B", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save(joinpath(@__DIR__, "figure_si_S6_1.png"), fig)


    fig4 = Figure(resolution=(1600*resfactor, 1000*resfactor))

    ax_plots = Axis(fig4[1:2, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_plots = rotr90(load(joinpath(@__DIR__, "screenshot_plots_cropped_annot.png")))
    hidedecorations!(ax_plots)
    image!(ax_plots, img_plots, aspect = DataAspect())

    #ax_plots_11 = Axis(fig4[3:5, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    #img_plots_11 = rotr90(load(joinpath(@__DIR__, "add_plots.png")))
    #hidedecorations!(ax_plots_11)
    #image!(ax_plots_11, img_plots_11, aspect = DataAspect())

    #ax_plots_12 = Axis(fig4[4, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    #img_plots_12 = rotr90(load(joinpath(@__DIR__, "screenshot_plots_oddsratio_dist_2.png")))
    #hidedecorations!(ax_plots_12)
    #image!(ax_plots_12, img_plots_12, aspect = DataAspect())

    #ax_plots_21 = Axis(fig4[4, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    #img_plots_21 = rotr90(load(joinpath(@__DIR__, "node_degree_plots.png")))
    #hidedecorations!(ax_plots_21)
    #image!(ax_plots_21, img_plots_21, aspect = DataAspect())

    #ax_plots_22 = Axis(fig4[5, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    #img_plots_22 = rotr90(load(joinpath(@__DIR__, "screenshot_plots_degree_dist_2.png")))
    #hidedecorations!(ax_plots_22)
    #image!(ax_plots_22, img_plots_22, aspect = DataAspect())

    #ax_plots_31 = Axis(fig4[5, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    #img_plots_31 = rotr90(load(joinpath(@__DIR__, "annotation_plots.png")))
    #hidedecorations!(ax_plots_31)
    #image!(ax_plots_31, img_plots_31, aspect = DataAspect())

    #ax_plots_32 = Axis(fig4[6, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    #img_plots_32 = rotr90(load(joinpath(@__DIR__, "screenshot_plots_annotations_2.png")))
    #hidedecorations!(ax_plots_32)
    #image!(ax_plots_32, img_plots_32, aspect = DataAspect())

    #Label(fig4[1,1, TopLeft()], "C", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    #Label(fig4[3,1, TopLeft()], "D", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    #Label(fig4[4,1, TopLeft()], "E", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    #Label(fig4[5,1, TopLeft()], "F", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save(joinpath(@__DIR__, "figure_si_S6_2_1.png"), fig4)


    fig5 = Figure(resolution=(1200*resfactor, 800*resfactor))

    ax_summary = Axis(fig5[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_summary = rotr90(load(joinpath(@__DIR__, "screenshot_summary_cropped_annot.png")))
    hidedecorations!(ax_summary)
    image!(ax_summary, img_summary, aspect = DataAspect())

    Label(fig5[1,1, TopLeft()], "G", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save(joinpath(@__DIR__, "figure_si_S6_3.png"), fig5)

end
si_plot()

using CairoMakie, FileIO

function si_plot()
    resfactor = 1.0
    fig = Figure(resolution=(1200*resfactor, 1200*resfactor))


    gb = fig[1,1] = GridLayout()

    ax_graph = Axis(gb[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(@__DIR__, "4A.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())

    ax_graph2 = Axis(gb[1, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph2 = rotr90(load(joinpath(@__DIR__, "4B.png")))
    hidedecorations!(ax_graph2)
    image!(ax_graph2, img_graph2, aspect = DataAspect())

    ga = fig[2,1] = GridLayout()

    ax_odds_fisher = Axis(ga[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_odds_fisher = rotr90(load(joinpath(@__DIR__, "oddsratio_fisher.png")))
    hidedecorations!(ax_odds_fisher)
    image!(ax_odds_fisher, img_odds_fisher, aspect = DataAspect())

    ax_odds_bp = Axis(ga[1, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_odds_fisher = rotr90(load(joinpath(@__DIR__, "oddsratio_bp.png")))
    hidedecorations!(ax_odds_bp)
    image!(ax_odds_bp, img_odds_fisher, aspect = DataAspect())

    Label(gb[1,1, TopLeft()], "A", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(gb[1,2, TopLeft()], "B", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,1, TopLeft()], "C", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,2, TopLeft()], "D", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save(joinpath(@__DIR__, "figure_si_S3.svg"), fig)

end
si_plot()

function odds_ratio_hists!(ax_fisher, ax_bp, interact::Interactions, plotting_fdr_level::Float64)
    valid_index = interact.edges.odds_ratio .>= 0.0
    logodds = log.(interact.edges.odds_ratio[valid_index])
    infindex = isfinite.(logodds)
    hist!(ax_fisher, logodds[infindex], bins=-8:0.1:15, label="all")
    sigfish_index = interact.edges.fisher_fdr[valid_index] .<= plotting_fdr_level
    hist!(ax_fisher, logodds[sigfish_index .& infindex], bins=-8:0.1:15, label="fdr <= $plotting_fdr_level")

    all_index = interact.edges.bp_fdr[valid_index] .<= 1.0
    hist!(ax_bp, logodds[all_index .& infindex], bins=-8:0.1:15, label="all")
    sigpred_index = interact.edges.bp_fdr[valid_index] .<= plotting_fdr_level
    hist!(ax_bp, logodds[sigpred_index .& infindex], bins=-8:0.1:15, label="fdr <= $plotting_fdr_level")
end

function plot_figure_s3(assets_folder::String, interact::Interactions)
    resfactor = 1.0
    fig = Figure(resolution=(1200*resfactor, 800*resfactor))

    gb = fig[1,1] = GridLayout()

    ax_graph = Axis(gb[1, 1], title="OppZ basepairing", ylabel="count", xlabel="position")#, leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    aggregation_plot!(ax_graph, interact, "OppZ", 0.3)
    axislegend(ax_graph)
    #hidexdecorations!(ax_graph, grid = false)
    #hidespines!(ax_graph, :l, :t, :r)

    ax_graph2 = Axis(gb[1, 2], title="VC2770 basepairing", ylabel="count", xlabel="position")#, leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    aggregation_plot!(ax_graph2, interact, "VC2770", 0.3)
    axislegend(ax_graph2)
    #hidexdecorations!(ax_graph2, grid = false)
    #hidespines!(ax_graph2, :l, :t, :r)

    ga = fig[2,1] = GridLayout()

    ax_odds_fisher = Axis(ga[1, 1], title="Odds ratio histogram for Fisher test")
    #    leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    #img_odds_fisher = rotr90(load(joinpath(assets_folder, "oddsratio_fisher.png")))
    #hidedecorations!(ax_odds_fisher)
    #image!(ax_odds_fisher, img_odds_fisher, aspect = DataAspect())

    ax_odds_bp = Axis(ga[1, 2], title="Odds ratio histogram for basepairing test")
    #    leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    #img_odds_fisher = rotr90(load(joinpath(assets_folder, "oddsratio_bp.png")))
    #hidedecorations!(ax_odds_bp)
    #image!(ax_odds_bp, img_odds_fisher, aspect = DataAspect())

    odds_ratio_hists!(ax_odds_fisher, ax_odds_bp, interact, 0.1)
    axislegend(ax_odds_fisher)
    axislegend(ax_odds_bp)

    Label(gb[1,1, TopLeft()], "A", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(gb[1,2, TopLeft()], "B", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,1, TopLeft()], "C", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,2, TopLeft()], "D", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_s3.svg", fig)

end

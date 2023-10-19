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

function plot_figure_s3(assets_folder::String, interact::Interactions, interact2::Interactions)
    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 1000*resfactor))

    ax_graph = Axis(fig[1, 1], title="OppZ basepairing", ylabel="count", xlabel="position")#, leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    aggregation_plot!(ax_graph, interact, "OppZ", 0.3)
    axislegend(ax_graph)

    ax_graph2 = Axis(fig[1, 2], title="VC2770 basepairing", ylabel="count", xlabel="position")#, leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    aggregation_plot!(ax_graph2, interact, "VC2770", 0.3)
    axislegend(ax_graph2)

    ax_pos = Axis(fig[2, 1], title="position in chimera", ylabel="count", xticks = (1:2, ["as RNA1", "as RNA2"]))
    ax_partner = Axis(fig[2, 2], title="partners of sRNAs", ylabel="average # of partners", xticks = (1:2, ["as RNA1", "as RNA2"]))

    srnas = unique(vcat(findall(interact.nodes.type .== "sRNA"), findall(interact2.nodes.type .== "sRNA")))
    is_sponge = falses(length(srnas))
    is_sponge2 = falses(length(srnas))
    counter_partners = zeros(Int, length(srnas))
    counter_partners2 = zeros(Int, length(srnas))
    for (i,srna) in enumerate(srnas)
        if srna <= length(interact.nodes.name)
            c1 = count_ligation_sites_as1(srna, interact, 30, 0.25)
            c2 = count_ligation_sites_as2(srna, interact, 30, 0.25)
            s1 = sum(v for pointvals in values(c1) for v in values(pointvals); init=0)
            s2 = sum(v for pointvals in values(c2) for v in values(pointvals); init=0)
            s1 > s2 && (is_sponge[i] = true)
            se = Set{String}()
            for k in keys.(values(c1))
                union!(se, Set(k))
            end
            for k in keys.(values(c2))
                union!(se, Set(k))
            end
            counter_partners[i] = length(se)
        end
        if srna <= length(interact2.nodes.name)
            c1 = count_ligation_sites_as1(srna, interact2, 30, 0.25)
            c2 = count_ligation_sites_as2(srna, interact2, 30, 0.25)
            s1 = sum(v for pointvals in values(c1) for v in values(pointvals); init=0)
            s2 = sum(v for pointvals in values(c2) for v in values(pointvals); init=0)
            s1 > s2 && (is_sponge2[i] = true)
            se = Set{String}()
            for k in keys.(values(c1))
                union!(se, Set(k))
            end
            for k in keys.(values(c2))
                union!(se, Set(k))
            end
            counter_partners2[i] = length(se)
        end
    end

    colors = Makie.wong_colors()
    mean_partners = [mean(counter_partners[is_sponge]), mean(counter_partners2[is_sponge2]), mean(counter_partners[.!is_sponge]), mean(counter_partners2[.!is_sponge2])]
    sd_partners = [std(counter_partners[is_sponge]), std(counter_partners2[is_sponge2]), std(counter_partners[.!is_sponge]), std(counter_partners2[.!is_sponge2])]

    exp = [1,1,2,2]
    groups = [1,2,1,2]

    barplot!(ax_partner, exp, mean_partners, dodge=groups, color=colors[groups])
    errorbars!(ax_partner, [0.8, 1.2, 1.8, 2.2], mean_partners, sd_partners, whiskerwidth=5)

    labels = ["HCD", "LCD"]
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    axislegend(ax_partner, elements, labels, position=:lt)

    count_sponge = [sum(is_sponge), sum(is_sponge2), sum(.!is_sponge), sum(.!is_sponge2)]

    barplot!(ax_pos, exp, count_sponge, dodge=groups, color=colors[groups])

    axislegend(ax_pos, elements, labels, position=:lt)

    ax_odds_fisher = Axis(fig[3, 1], title="Odds ratio histogram for Fisher test", ylabel="count", xlabel="log(odds ratio)")
    ax_odds_bp = Axis(fig[3, 2], title="Odds ratio histogram for basepairing test", ylabel="count", xlabel="log(odds ratio)")

    odds_ratio_hists!(ax_odds_fisher, ax_odds_bp, interact, 0.1)

    axislegend(ax_odds_fisher)
    axislegend(ax_odds_bp)

    Label(fig[1,1, TopLeft()], "a", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[1,2, TopLeft()], "b", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,1, TopLeft()], "c", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,2, TopLeft()], "d", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,1, TopLeft()], "e", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,2, TopLeft()], "f", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_s3.svg", fig)
    save("figure_s3.png", fig, px_per_unit = 2)
end

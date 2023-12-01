function odds_ratio_hists!(ax_fisher, ax_bp, interact::InteractionsNew, plotting_fdr_level::Float64)
    colors = ((:tomato4, 0.6), (:forestgreen, 0.7))
    valid_index = interact.edges.odds_ratio .>= 0.0
    logodds = log.(interact.edges.odds_ratio[valid_index])
    infindex = isfinite.(logodds)
    hist!(ax_fisher, logodds[infindex], bins=-8:0.15:15, label="all", color=colors[1])
    sigfish_index = interact.edges.fisher_fdr[valid_index] .<= plotting_fdr_level
    hist!(ax_fisher, logodds[sigfish_index .& infindex], bins=-8:0.15:15, label="Fisher FDR <= $plotting_fdr_level", color=colors[2])

    all_index = interact.edges.bp_fdr[valid_index] .<= 1.0
    hist!(ax_bp, logodds[all_index .& infindex], bins=-8:0.15:15, label="all", color=colors[1])
    sigpred_index = interact.edges.bp_fdr[valid_index] .<= plotting_fdr_level
    hist!(ax_bp, logodds[sigpred_index .& infindex], bins=-8:0.15:15, label="basepairing FDR <= $plotting_fdr_level", color=colors[2])
end

function extract_counts(interacts, srna_type, mrna_type, fdr_cut)
    counts =[
        begin
            srna1_index = interact.nodes.type[interact.edges.src] .== srna_type
            srna2_index = interact.nodes.type[interact.edges.dst] .== srna_type
            mrna1_index = interact.nodes.type[interact.edges.src] .== mrna_type
            mrna2_index = interact.nodes.type[interact.edges.dst] .== mrna_type
            srna1_index = trues(length(interact.edges.dst))
            srna2_index = trues(length(interact.edges.dst))
            mrna1_index = trues(length(interact.edges.dst))
            mrna2_index = trues(length(interact.edges.dst))
            combined_index = (srna1_index .& mrna2_index) .| (mrna1_index .& srna2_index)
            pairs = Dict{Tuple{String, String}, Bool}()
            pairs_sig = Dict{Tuple{String, String}, Bool}()
            for row in eachrow(@view interact.edges[srna1_index .& mrna2_index, :])
                pair = (interact.nodes.name[row.src],interact.nodes.name[row.dst])
                pairs[pair]=!isnan(row.bp_fdr)
                pairs_sig[pair]=row.bp_fdr <= fdr_cut
            end
            for row in eachrow(@view interact.edges[mrna1_index .& srna2_index, :])
                pair = (interact.nodes.name[row.src],interact.nodes.name[row.dst])
                pairs[pair]=!isnan(row.bp_fdr)
                pairs_sig[pair]=row.bp_fdr <= fdr_cut
            end
            unique_pairs = Dict{Tuple{String, String}, Bool}()
            unique_pairs_sig = Dict{Tuple{String, String}, Bool}()
            for (pair, has_ligation) in pairs
                is_significant = pairs_sig[pair]
                revpair = (last(pair), first(pair))
                if revpair in keys(unique_pairs)
                    unique_pairs[revpair] = unique_pairs[revpair] | has_ligation
                    unique_pairs_sig[revpair] = unique_pairs_sig[revpair] | is_significant
                else
                    unique_pairs[pair] = has_ligation
                    unique_pairs_sig[pair] = is_significant
                end
            end
            [sum(interact.edges.nb_ints[combined_index]), length(pairs),
                length(unique_pairs), sum(values(unique_pairs)), sum(values(unique_pairs_sig))]
        end for interact in interacts
    ]
    return counts
end

function plot_figure_s4(interacts_vibrio, interacts_ecoli_rilseq, interacts_ecoli_clash,
        interacts_subtilis_ligrseq, interacts_salmonella, interacts_pseudomonas, interacts_epec, reporteds, fdr_cut)

    c_vib = extract_counts(interacts_vibrio, "sRNA", "CDS_UTRS", fdr_cut)
    c_eco_rilseq = extract_counts(interacts_ecoli_rilseq, "ncRNA", "CDS_UTR", fdr_cut)
    c_eco_clash = extract_counts(interacts_ecoli_clash, "ncRNA", "CDS_UTRS", fdr_cut)
    c_bsub_ligr = extract_counts(interacts_subtilis_ligrseq, "ncRNA", "CDS_UTRS", fdr_cut)
    c_salmonella = extract_counts(interacts_salmonella, "ncRNA", "CDS_UTRS", fdr_cut)
    c_pseudomonas = extract_counts(interacts_pseudomonas, "ncRNA", "CDS_UTRS", fdr_cut)
    c_epec = extract_counts(interacts_epec, "ncRNA", "CDS_UTRS", fdr_cut)

    norms = [10^4, 10^4, 10^3, 10^4, 10^4, 10^4, 10^4]

    colors = Makie.wong_colors()
    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 1400*resfactor))

    ax_odds_fisher = Axis(fig[1, 1], title="Odds ratio histogram for Fisher test", ylabel="count / 10^2", xlabel="log(odds ratio)",
        yticks=(500:500:1500, ["5", "10", "15"]))
    ax_odds_bp = Axis(fig[1, 2], title="Odds ratio histogram for complementarity test", ylabel="count / 10^2", xlabel="log(odds ratio)",
        yticks=(200:200:800, ["2", "4", "6", "8"]))

    odds_ratio_hists!(ax_odds_fisher, ax_odds_bp, interacts_vibrio[1], 0.1)

    axislegend(ax_odds_fisher, position=:rt)
    axislegend(ax_odds_bp, position=:rt)

    ax_vib = Axis(fig[2,1], title=rich("RIL-seq ", rich("V. cholerae"; font=:bold_italic)),
        xlabel="dataset", ylabel="count / 10^4", xticks=(1:2, ["LCD", "HCD"]))
    ax_eco_rilseq = Axis(fig[3,1], title=rich("RIL-seq ", rich("E. coli"; font=:bold_italic)),
        xlabel="dataset", ylabel="count / 10^4", xticks=(1:3, ["exponential", "stationary", "iron limit"]))
    ax_epec = Axis(fig[3,2], title=rich("RIL-seq ", rich("EPEC"; font=:bold_italic)),
        xlabel="dataset", ylabel="count / 10^4", xticks=(1:2, ["wild-type", "T3SS_BSP"]))
    ax_salmonella = Axis(fig[4,1], title=rich("RIL-seq ", rich("S. enterica"; font=:bold_italic)),
        xlabel="dataset", ylabel="count / 10^4", xticks=(1:1, ["wild-type"]))
    ax_pseudomonas = Axis(fig[4,2], title=rich("RIL-seq ", rich("P. aeruginosa"; font=:bold_italic)),
        xlabel="dataset", ylabel="count / 10^4", xticks=(1:2, ["exponential", "stationary"]))
    ax_eco_clash = Axis(fig[5,1], title=rich("CLASH ", rich("E. coli"; font=:bold_italic)),
        xlabel="dataset", ylabel="count / 10^3", xticks=(1:3, ["exponential", "transition", "stationary"]), yticks=1:5)
    ax_bsub_ligr = Axis(fig[5,2], title=rich("LIGR-seq ", rich("B. subtilis"; font=:bold_italic)),
        xlabel="dataset", ylabel="count / 10^4", xticks=(1:1, ["combined"]))


    for (counts, ax, reported, n) in zip((c_vib, c_eco_rilseq, c_eco_clash, c_bsub_ligr, c_salmonella, c_pseudomonas, c_epec),
    					(ax_vib, ax_eco_rilseq, ax_eco_clash, ax_bsub_ligr, ax_salmonella, ax_pseudomonas, ax_epec),
            reporteds, norms)
        c = vcat([[r, count[2:5]...] ./ n for (r,count) in zip(reported, counts)]...)
        exp = vcat([[i, i, i, i, i] for i in 1:length(counts)]...)
        groups = vcat([[1, 2, 3, 4, 5] for i in 1:length(counts)]...)
        barplot!(ax, exp, c, dodge=groups, color=colors[groups])
        hidexdecorations!(ax, label = false, ticklabels = false, ticks = false,
            grid = true, minorgrid = false, minorticks = false)
    end
    labels = [["reported"], ["total", "unique", "with ligation", "significant"]]
    elements = [[PolyElement(polycolor = colors[1])], [PolyElement(polycolor = colors[i]) for i in 2:5]]

    l = Legend(fig[2,2], elements, labels, ["original study", "ChimericFragments"])
    l.orientation = :horizontal
    l.nbanks = 2
    l.tellwidth = false
    l.tellheight = false
    l.groupgap = 50

    Label(fig[1,1, TopLeft()], "a", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[1,2, TopLeft()], "b", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,1, TopLeft()], "c", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[5,2, TopLeft()], "i", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,1, TopLeft()], "d", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,2, TopLeft()], "e", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[4,1, TopLeft()], "f", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[4,2, TopLeft()], "g", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[5,1, TopLeft()], "h", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_E6.svg", fig)
    save("figure_E6.png", fig, px_per_unit = 2)
end

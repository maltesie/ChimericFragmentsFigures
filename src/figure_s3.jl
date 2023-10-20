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

function hdi_index(counts::Vector{Int}; probmasscut=0.8)
    probMassVec = counts/sum(counts)
    sortedProbMass = sort(probMassVec; rev=true)
    HDIheightIdx = findfirst(cumsum(sortedProbMass) .>= probmasscut)
    HDIheight = sortedProbMass[HDIheightIdx]
    probMassVec .>= HDIheight
end

function count_intervals(index::BitVector; connect_within=1)
    my_index = copy(index)
    for i in (connect_within+1):(length(index)-connect_within)
        my_index[i] = my_index[i] | (any(my_index[(i-connect_within):(i-1)]) & any(my_index[(i+1):(i+connect_within)]))
    end
    sum(diff(my_index) .== 1)
end

function values_of_counts(counts::Dict{Int, Dict{String, Int}})
    ma, mi = maximum(keys(counts)), minimum(keys(counts))
    vals = zeros(Int, ma-mi+1)
    for (i, v) in enumerate(mi:ma)
        v in keys(counts) && (vals[i] = sum(values(counts[v]); init=0))
    end
    vals
end

countofnints(vals::Vector{Int}, n::Int) = [sum(vals .== i) for i in 1:n]

function plot_figure_s3(assets_folder::String, interact::Interactions, interact2::Interactions;
        bp_interval_len=30, fdr_cut=0.25, probmass_min=0.3, probmass_max=0.7, min_partners=3, min_ligpoints=10, connect_within=1)
    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 1200*resfactor))

    ax_graph = Axis(fig[1, 1], title="OppZ basepairing", ylabel="count", xlabel="position")
    aggregation_plot!(ax_graph, interact, "OppZ", fdr_cut)
    axislegend(ax_graph)

    ax_graph2 = Axis(fig[1, 2], title="oppB basepairing", ylabel="count", xlabel="position")
    aggregation_plot!(ax_graph2, interact, "oppB", fdr_cut)
    axislegend(ax_graph2)

    ax_pos = Axis(fig[3, 1], title="position in chimera", ylabel="count", xticks = (1:2, ["RNA1", "RNA2"]))
    ax_partner = Axis(fig[3, 2], title="partners of sRNAs", ylabel="# of partners", xticks = (1:2, ["RNA1", "RNA2"]))

    ax_regions_s = Axis(fig[2, 1], title="regions in RNA1 sRNAs", ylabel="count", xlabel="# of regions", xticks = 1:100)
    ax_regions = Axis(fig[2, 2], title="regions in RNA2 sRNAs", ylabel="count", xlabel="# of regions", xticks = 1:100)

    srnas = unique(vcat(findall(interact.nodes.type .== "sRNA"), findall(interact2.nodes.type .== "sRNA")))
    is_sponge = falses(length(srnas))
    is_sponge2 = falses(length(srnas))
    counter_partners = zeros(Int, length(srnas))
    counter_partners2 = zeros(Int, length(srnas))
    counter_reads = zeros(Int, (length(srnas),2))
    counter_reads2 = zeros(Int, (length(srnas),2))
    counter_regions = zeros(Int, (length(srnas),2))
    counter_regions2 = zeros(Int, (length(srnas),2))

    for (i,srna) in enumerate(srnas)
        if srna <= length(interact.nodes.name)
            c1 = count_ligation_sites_as1(srna, interact, bp_interval_len, fdr_cut)
            c2 = count_ligation_sites_as2(srna, interact, bp_interval_len, fdr_cut)
            #s1 = sum(v for pointvals in values(c1) for v in values(pointvals); init=0)
            #s2 = sum(v for pointvals in values(c2) for v in values(pointvals); init=0)
            #s1 > s2 && (is_sponge[i] = true)
            se1 = Set{String}()
            for k in keys.(values(c1))
                union!(se1, Set(k))
            end
            se2 = Set{String}()
            for k in keys.(values(c2))
                union!(se2, Set(k))
            end

            length(se1) > length(se2) && (is_sponge[i] = true)

            counter_partners[i] = length(union(se1, se2))

            if length(c1) > 0
                vals1 = values_of_counts(c1)
                counter_reads[i, 1] = maximum(vals1)
                n_intervals1 = maximum(count_intervals(hdi_index(vals1; probmasscut=pc); connect_within=connect_within) for pc in probmass_min:0.1:probmass_max)
                counter_regions[i, 1] = n_intervals1
            end
            if length(c2) > 0
                vals2 = values_of_counts(c2)
                counter_reads[i, 2] = maximum(vals2)
                n_intervals2 = maximum(count_intervals(hdi_index(vals2; probmasscut=pc); connect_within=connect_within) for pc in probmass_min:0.1:probmass_max)
                counter_regions[i, 2] = n_intervals2
            end
        end

        if srna <= length(interact2.nodes.name)
            c1 = count_ligation_sites_as1(srna, interact2, bp_interval_len, fdr_cut)
            c2 = count_ligation_sites_as2(srna, interact2, bp_interval_len, fdr_cut)
            #s1 = sum(v for pointvals in values(c1) for v in values(pointvals); init=0)
            #s2 = sum(v for pointvals in values(c2) for v in values(pointvals); init=0)
            #s1 > s2 && (is_sponge2[i] = true)
            se1 = Set{String}()
            for k in keys.(values(c1))
                union!(se1, Set(k))
            end
            se2 = Set{String}()
            for k in keys.(values(c2))
                union!(se2, Set(k))
            end
            length(se1) > length(se2) && (is_sponge2[i] = true)
            counter_partners2[i] = length(union(se1, se2))

            if length(c1) > 0
                vals1 = values_of_counts(c1)
                counter_reads2[i, 1] = maximum(vals1)
                n_intervals1 = maximum(count_intervals(hdi_index(vals1; probmasscut=pc); connect_within=connect_within) for pc in probmass_min:0.1:probmass_max)
                counter_regions2[i, 1] = n_intervals1
            end
            if length(c2) > 0
                vals2 = values_of_counts(c2)
                counter_reads2[i, 2] = maximum(vals2)
                n_intervals2 = maximum(count_intervals(hdi_index(vals2; probmasscut=pc); connect_within=connect_within) for pc in probmass_min:0.1:probmass_max)
                counter_regions2[i, 2] = n_intervals2
            end
        end
    end

    filterindex = vec((maximum(counter_reads; dims=2) .>= min_ligpoints) .& (counter_partners .>= min_partners))
    is_sponge = is_sponge[filterindex]
    counter_regions = counter_regions[filterindex, :]
    counter_partners = counter_partners[filterindex]

    filterindex2 = vec((maximum(counter_reads2; dims=2) .>= min_ligpoints) .& (counter_partners2 .>= min_partners))
    is_sponge2 = is_sponge2[filterindex2]
    counter_regions2 = counter_regions2[filterindex2, :]
    counter_partners2 = counter_partners2[filterindex2]

    colors = Makie.wong_colors()
    mean_partners = [mean(counter_partners[is_sponge]), mean(counter_partners2[is_sponge2]), mean(counter_partners[.!is_sponge]), mean(counter_partners2[.!is_sponge2])]
    sd_partners = [std(counter_partners[is_sponge]), std(counter_partners2[is_sponge2]), std(counter_partners[.!is_sponge]), std(counter_partners2[.!is_sponge2])]

    len_p1_s = length(counter_partners[is_sponge])
    len_p2_s = length(counter_partners2[is_sponge2])
    len_p1 = length(counter_partners[.!is_sponge])
    len_p2 = length(counter_partners2[.!is_sponge2])
    ordered_partners = vcat(counter_partners[is_sponge], counter_partners2[is_sponge2], counter_partners[.!is_sponge], counter_partners2[.!is_sponge2])
    ordered_exp = vcat(repeat([1], len_p1_s), repeat([1], len_p2_s), repeat([2], len_p1), repeat([2], len_p2))
    ordered_groups = vcat(repeat([1], len_p1_s), repeat([2], len_p2_s), repeat([1], len_p1), repeat([2], len_p2))

    #barplot!(ax_partner, exp, mean_partners, dodge=groups, color=colors[groups])
    #errorbars!(ax_partner, [0.8, 1.2, 1.8, 2.2], mean_partners, sd_partners, whiskerwidth=5)
    #violin!(ax_partner, ordered_exp, ordered_partners, dodge=ordered_groups, color=colors[ordered_groups], datalimits=(0,500))
    boxplot!(ax_partner, ordered_exp, ordered_partners, dodge=ordered_groups, color=colors[ordered_groups], show_outliers=false)

    labels = ["HCD", "LCD"]
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    axislegend(ax_partner, elements, labels, position=:lt)

    count_sponge = [sum(is_sponge), sum(is_sponge2), sum(.!is_sponge), sum(.!is_sponge2)]
    exp = [1,1,2,2]
    groups = [1,2,1,2]

    barplot!(ax_pos, exp, count_sponge, dodge=groups, color=colors[groups])

    axislegend(ax_pos, elements, labels, position=:lt)

    n_s = max(maximum(counter_regions[is_sponge, 1]), maximum(counter_regions2[is_sponge2, 1]))

    #count_regions = vcat(counter_regions[is_sponge, 1], counter_regions2[is_sponge2, 1], counter_regions[.!is_sponge, 2], counter_regions2[.!is_sponge2, 2])
    r1_s = countofnints(counter_regions[is_sponge, 1], n_s)
    r2_s = countofnints(counter_regions2[is_sponge2, 1], n_s)

    count_regions_s= vcat(r1_s, r2_s)
    #regions_exp_s = vcat(repeat([1], 10), repeat([1], len_p2_s), repeat([2], len_p1), repeat([2], len_p2))
    regions_groups_s = vcat(repeat([1], n_s), repeat([2], n_s))

    barplot!(ax_regions_s, vcat(1:n_s, 1:n_s), count_regions_s, dodge=regions_groups_s, color=colors[regions_groups_s])
    axislegend(ax_regions_s, elements, labels, position=:rt)

    n = max(maximum(counter_regions[.!is_sponge, 2]), maximum(counter_regions2[.!is_sponge2, 2]))

    r1 = countofnints(counter_regions[.!is_sponge, 2], n)
    r2 = countofnints(counter_regions2[.!is_sponge2, 2], n)

    count_regions= vcat(r1, r2)
    #regions_exp_s = vcat(repeat([1], 10), repeat([1], len_p2_s), repeat([2], len_p1), repeat([2], len_p2))
    regions_groups = vcat(repeat([1], n), repeat([2], n))

    barplot!(ax_regions, vcat(1:n, 1:n), count_regions, dodge=regions_groups, color=colors[regions_groups])
    axislegend(ax_regions, elements, labels, position=:rt)

    ax_odds_fisher = Axis(fig[4, 1], title="Odds ratio histogram for Fisher test", ylabel="count", xlabel="log(odds ratio)")
    ax_odds_bp = Axis(fig[4, 2], title="Odds ratio histogram for basepairing test", ylabel="count", xlabel="log(odds ratio)")

    odds_ratio_hists!(ax_odds_fisher, ax_odds_bp, interact, 0.1)

    axislegend(ax_odds_fisher)
    axislegend(ax_odds_bp)

    Label(fig[1,1, TopLeft()], "a", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[1,2, TopLeft()], "b", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,1, TopLeft()], "c", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,2, TopLeft()], "d", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,1, TopLeft()], "e", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,2, TopLeft()], "f", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[4,1, TopLeft()], "g", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[4,2, TopLeft()], "h", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_s3.svg", fig)
    save("figure_s3.png", fig, px_per_unit = 2)
end

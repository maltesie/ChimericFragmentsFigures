function count_ligation_sites_as1(node_id::Int, interact::InteractionsNew, bp_len::Int, max_fdr::Float64; fisher_fdr=1.0)
    df = interact.edges
    counts = Dict{Int, Dict{String,Int}}()
    isnegative = interact.nodes.strand[node_id] == '-'
    for row in eachrow(@view df[df.src .== node_id, [:src, :dst, :fisher_fdr]])
        partner = row.dst
        row.fisher_fdr > fisher_fdr && continue
        ligation_points = interact.edgestats[(node_id, partner)][3]
        length(ligation_points) > 0 || continue
        n = interact.nodes.name[partner]
        fdr = adjust(PValues([interact.bpstats[(row.src, first(p), row.dst, last(p))][1] for p in keys(ligation_points)]), BenjaminiHochberg())
        for ((p, c), f) in zip(ligation_points, fdr)
            f < max_fdr || continue
            for pl in interact.bpstats[(row.src, first(p), row.dst, last(p))][2]:interact.bpstats[(row.src, first(p), row.dst, last(p))][3]
                p1 = p[1] + (isnegative ? 1 : -1) * (bp_len - pl)
                if p1 in keys(counts)
                    n in keys(counts[p1]) ? counts[p1][n] += c : counts[p1][n] = c
                else
                    counts[p1] = Dict(n=>c)
                end
            end
        end
    end
    return counts
end

function count_ligation_sites_as2(node_id::Int, interact::InteractionsNew, bp_len::Int, max_fdr::Float64; fisher_fdr=1.0)
    df = interact.edges
    counts = Dict{Int, Dict{String,Int}}()
    isnegative = interact.nodes.strand[node_id] == '-'
    for row in eachrow(@view df[df.dst .== node_id, [:src, :dst, :fisher_fdr]])
        partner = row.src
        row.fisher_fdr > fisher_fdr && continue
        ligation_points = interact.edgestats[(partner, node_id)][3]
        length(ligation_points) > 0 || continue
        n = interact.nodes.name[partner]
        fdr = adjust(PValues([interact.bpstats[(row.src, first(p), row.dst, last(p))][1] for p in keys(ligation_points)]), BenjaminiHochberg())
        for ((p, c), f) in zip(ligation_points, fdr)
            f < max_fdr || continue
            for pl in interact.bpstats[(row.src, first(p), row.dst, last(p))][4]:interact.bpstats[(row.src, first(p), row.dst, last(p))][5]
                p2 = p[2] + (isnegative ? -1 : 1) * (bp_len - pl) - 1
                if p2 in keys(counts)
                    n in keys(counts[p2]) ? counts[p2][n] += c : counts[p2][n] = c
                else
                    counts[p2] = Dict(n=>c)
                end
            end
        end
    end
    return counts
end

function aggregation_plot!(ax::Axis, interact::InteractionsNew, n::String, max_fdr::Float64; norm=1.0)
    idx = findfirst(interact.nodes.name .== n)
    cds = interact.nodes.cds[idx] == 0 ? (interact.nodes.strand[idx] == '-' ? interact.nodes.right[idx] : interact.nodes.left[idx]) : interact.nodes.cds[idx]
    st = interact.nodes.strand[idx]
    count1 = count_ligation_sites_as1(idx, interact, 30, max_fdr)
    count2 = count_ligation_sites_as2(idx, interact, 30, max_fdr)
    (length(count1) == 0) && (length(count2) == 0) && return nothing
    for (color, label, count) in zip((RGBAf(0.388, 0.431, 0.98, 0.6), RGBAf(0.937, 0.333, 0.231, 0.6)), ("RNA1", "RNA2"), (count1, count2))
        ligationpoints = count
        isempty(ligationpoints) && continue
        kv = collect(sort(ligationpoints, by=x->x[1]))
        positions, counts = first.(kv), sum.(values.(last.(kv)))
        allpositions = length(positions) > 0 ? collect(minimum(positions):maximum(positions)) : Int[]
        pindex = in.(allpositions, Ref(positions))
        indextrans = [pindex[i] ? sum(view(pindex, 1:i)) : 0 for i in eachindex(allpositions)]
        allcounts = [pindex[i] ? counts[indextrans[i]] : 0 for i in eachindex(allpositions)]
        allcounts_normed = allcounts ./ norm
        if st == '+'
            allpositions .-= cds
        else
            allpositions = reverse(-1.0 .* (allpositions .- cds))
            reverse!(allcounts_normed)
        end
        allpositions .+= 1
        band!(ax, allpositions, zeros(length(allcounts)), allcounts_normed, label=label, color=color)
    end
end

function seed_correlation_plot!(ax::Axis, interact::InteractionsNew, n::String, range1::UnitRange, range2::UnitRange, hn1::String, hn2::String, max_fdr::Float64)
    idx = findfirst(interact.nodes.name .== n)
    cds = interact.nodes.cds[idx] == 0 ? (interact.nodes.strand[idx] == '-' ? interact.nodes.right[idx] : interact.nodes.left[idx]) : interact.nodes.cds[idx]
    st = interact.nodes.strand[idx]
    count1 = count_ligation_sites_as1(idx, interact, 30, max_fdr)
    count2 = count_ligation_sites_as2(idx, interact, 30, max_fdr)
    (length(count1) == 0) && (length(count2) == 0) && return nothing
    r1 = st == '+' ? range1 .+ (cds - 1) : (cds + 1) .- range1
    r2 = st == '+' ? range2 .+ (cds - 1) : (cds + 1) .- range2
    for (color, label, count) in zip((RGBAf(0.388, 0.431, 0.98, 0.7), RGBAf(0.937, 0.333, 0.231, 0.7)), ("RNA1", "RNA2"), (count1, count2))
        ligationpoints = count
        isempty(ligationpoints) && continue

        region1 = Dict(nid=>0 for p in values(count) for nid in keys(p))
        region2 = Dict(nid=>0 for p in values(count) for nid in keys(p))

        for nid in keys(region1)
            for (i,p) in enumerate(r1)
                if (p in keys(count)) && (nid in keys(count[p])) && (count[p][nid] > region1[nid])
                    region1[nid] = count[p][nid]
                end
            end
            for (i,p) in enumerate(r2)
                if (p in keys(count)) && (nid in keys(count[p])) && (count[p][nid] > region2[nid])
                    region2[nid] = count[p][nid]
                end
            end
        end
        ks = collect(keys(region1))
        scatter!(ax, log.([region1[k] for k in ks] .+ 1), log.([region2[k] for k in ks] .+ 1), label=label, color=color)
        #println(label, ": ", corspearman([region1[k] for k in ks], [region2[k] for k in ks]))
        if label == "RNA2"
            if hn1 in keys(region1)
                text!(ax, log(region1[hn1]+1), log(region2[hn1]+1); text=hn1)
            end
            if hn2 in keys(region1)
                text!(ax, log(region1[hn2]+1), log(region2[hn2]+1); text=hn2)
            end
        end
    end
end

function plot_figure_5(assets_folder::String, interact::InteractionsNew)

    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 800*resfactor))

    Label(fig[1,1, TopLeft()], "a", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(fig[1,2, TopLeft()], "b", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    Label(fig[2,1, TopLeft()], "c", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,2, TopLeft()], "d", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    colors = reverse(ColorSchemes.tofino10.colors)[[4,3,2,1]]
    #colors = ColorSchemes.berlin25.colors[[7,5,3,1]]
    #colors = ColorSchemes.berlin25.colors[[3,5,7,9]]
    #colors = Makie.wong_colors()
    #colors = ColorSchemes.grayC10.colors[[2,4,6,8]]

    ax3 = Axis(fig[1, 1], title="FarS basepairing", ylabel="count", xlabel="position")
    aggregation_plot!(ax3, interact, "FarS", 0.25)
    band!(ax3, [30., 40.], [-150., -150.], [0., 0.], color=colors[3], label="seed 1")
    band!(ax3, [50., 60.], [-150., -150.], [0., 0.], color=colors[4], label="seed 2")
    ylims!(ax3, (-150, 2500))
    axislegend(ax3, position=:rt)


    ax_seed_corr = Axis(fig[1, 2], title="FarS seeds correlation", ylabel="log count in seed 2", xlabel="log count in seed 1", xticks=2:2:6)
    seed_correlation_plot!(ax_seed_corr, interact, "FarS", 30:39, 44:50, "VC1043", "VCA0848", 0.25)
    xlims!(ax_seed_corr, (-0.5, 7))
    ylims!(ax_seed_corr, (-0.5, 5))
    Legend(fig[1,3], ax_seed_corr)

    img4 = rotr90(load(joinpath(assets_folder, "basepairing.png")))
    ax4 = Axis(fig[2, 1], title="basepairing predictions")#, aspect = DataAspect())
    hidedecorations!(ax4)
    image!(ax4, img4, aspect = DataAspect())

    df_plate = DataFrame(CSV.File(joinpath(assets_folder, "platereader2.csv")))

    ax5 = Axis(fig[2,2], title = "FarS reporter assay", ylabel="relative fluorescence [AU]",
        xticks = (1:2, [rich("vc1043"; font=:italic), rich("vca0848"; font=:italic)]))

    exp = [1,1,1,1,2,2,2,2]
    groups = [1,2,3,4,1,2,3,4]

    ctrl = Matrix(df_plate[:, [:VCA0848_1, :VCA0848_2, :VCA0848_3]])
    ctrl ./= mean(ctrl[1, :])

    spot = Matrix(df_plate[:, [:VC1043_1, :VC1043_2, :VC1043_3]])
    spot ./= mean(spot[1, :])

    mean_ctrl = vec(mean(ctrl; dims=2))
    sd_ctrl = vec(std(ctrl; dims=2))

    mean_spot = vec(mean(spot; dims=2))
    sd_spot = vec(std(spot; dims=2))

    barplot!(ax5, exp, vcat(mean_spot, mean_ctrl), dodge=groups, color=colors[groups])
    ctrlpos = [1.7, 1.9, 2.1, 2.3]
    scatter!(ax5, ctrlpos .- 0.05, vec(ctrl[:, 1]), color="black", markersize=5)
    scatter!(ax5, ctrlpos, vec(ctrl[:, 2]), color="black", markersize=5)
    scatter!(ax5, ctrlpos .+ 0.05, vec(ctrl[:, 3]), color="black", markersize=5)
    errorbars!(ax5, ctrlpos, mean_ctrl, sd_ctrl, whiskerwidth=5)
    spotpos = [0.7, 0.9, 1.1, 1.3]
    scatter!(ax5, spotpos .- 0.05, vec(spot[:, 1]), color="black", markersize=5)
    scatter!(ax5, spotpos, vec(spot[:, 2]), color="black", markersize=5)
    scatter!(ax5, spotpos .+ 0.05, vec(spot[:, 3]), color="black", markersize=5)
    errorbars!(ax5, spotpos, mean_spot, sd_spot, whiskerwidth=5)
    labels = df_plate.name
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    #axislegend(ax5, elements, labels, position=:rt)
    Legend(fig[2,3], elements, labels)
    hidexdecorations!(ax5, label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)

    ps = pvalue.(EqualVarianceTTest.([ctrl[1, :] for i in 2:size(ctrl)[1]], [ctrl[i, :] for i in 2:size(ctrl)[1]]))
    for (i, p) in enumerate(ps)
        i+=1
        if p > 0.05
            text!(ax5, "n.s.", position = (ctrlpos[i], mean_ctrl[i]+sd_ctrl[i] + 0.1), align = (:center, :center), fontsize=17)
        elseif p > 0.01
            text!(ax5, "⭑", position = (ctrlpos[i], mean_ctrl[i]+sd_ctrl[i] + 0.1), align = (:center, :center), fontsize=17)
        elseif p > 0.001
            text!(ax5, "⭑⭑", position = (ctrlpos[i], mean_ctrl[i]+sd_ctrl[i] + 0.1), align = (:center, :center), fontsize=17)
        else
            text!(ax5, "⭑⭑⭑", position = (ctrlpos[i], mean_ctrl[i]+sd_ctrl[i] + 0.1), align = (:center, :center), fontsize=17)
        end
    end

    ps = pvalue.(EqualVarianceTTest.([spot[1, :] for i in 2:size(spot)[1]], [spot[i, :] for i in 2:size(spot)[1]]))
    for (i, p) in enumerate(ps)
        i+=1
        if p > 0.05
            text!(ax5, "n.s.", position = (spotpos[i], mean_spot[i]+sd_spot[i] + 0.1), align = (:center, :center), fontsize=17)
        elseif p > 0.01
            text!(ax5, "⭑", position = (spotpos[i], mean_spot[i]+sd_spot[i] + 0.1), align = (:center, :center), fontsize=17)
        elseif p > 0.001
            text!(ax5, "⭑⭑", position = (spotpos[i], mean_spot[i]+sd_spot[i] + 0.1), align = (:center, :center), fontsize=17)
        else
            text!(ax5, "⭑⭑⭑", position = (spotpos[i], mean_spot[i]+sd_spot[i] + 0.1), align = (:center, :center), fontsize=17)
        end
    end

    #titlelayout = GridLayout(fig[0, 1], halign = :left, tellwidth = false)
    #Label(titlelayout[1, 1], "Fig. 5", halign = :left, fontsize=30)

    save("figure_5.pdf", fig)
    save("figure_5.png", fig, px_per_unit = 2)

end

function all_srna_aggregation_plots(interact_lcd::InteractionsNew, interact_hcd::InteractionsNew)

    resfactor = 1.
    fig = Figure(resolution=(800*resfactor, 1200*resfactor))
    srnas = sort(unique(vcat(interact_lcd.nodes.name[interact_lcd.nodes.type .== "sRNA"], interact_hcd.nodes.name[interact_hcd.nodes.type .== "sRNA"])))
    for (i,srna) in enumerate(srnas)
        idx = (i % 5 == 0) ? 5 : (i % 5)
        if srna in interact_lcd.nodes.name
            ax1 = Axis(fig[idx, 1], title="$srna LCD")
            aggregation_plot!(ax1, interact_lcd, srna, 0.25)
        end
        if srna in interact_hcd.nodes.name
            ax2 = Axis(fig[idx, 2], title="$srna HCD")
            aggregation_plot!(ax2, interact_hcd, srna, 0.25)
        end
        if i % 5 == 0
            save("figure_srnas_$(Int(floor(i/5))).pdf", fig, px_per_unit = 2)
            fig = Figure(resolution=(800*resfactor, 1200*resfactor))
            srnas = sort(unique(vcat(interact_lcd.nodes.name[interact_lcd.nodes.type .== "sRNA"], interact_hcd.nodes.name[interact_hcd.nodes.type .== "sRNA"])))
        end
    end
    #save("figure_srnas_99.svg", fig)
    save("figure_srnas_99.pdf", fig, px_per_unit = 2)
end
function count_ligation_sites_as1(node_id::Int, interact::Interactions, bp_len::Int, max_fdr::Float64)
    df = interact.edges
    counts = Dict{Int, Dict{String,Int}}()
    isnegative = interact.nodes.strand[node_id] == '-'
    for partner in df[df.src .== node_id, :dst]
        ligation_points = interact.edgestats[(node_id, partner)][3]
        length(ligation_points) > 0 || continue
        n = interact.nodes.name[partner]
        fdr = adjust(PValues([interact.bpstats[p][1] for p in keys(ligation_points)]), BenjaminiHochberg())
        for ((p, c), f) in zip(ligation_points, fdr)
            f < max_fdr || continue
            for pl in interact.bpstats[p][2]:interact.bpstats[p][3]
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

function count_ligation_sites_as2(node_id::Int, interact::Interactions, bp_len::Int, max_fdr::Float64)
    df = interact.edges
    counts = Dict{Int, Dict{String,Int}}()
    isnegative = interact.nodes.strand[node_id] == '-'
    for partner in df[df.dst .== node_id, :src]
        ligation_points = interact.edgestats[(partner, node_id)][3]
        length(ligation_points) > 0 || continue
        n = interact.nodes.name[partner]
        fdr = adjust(PValues([interact.bpstats[p][1] for p in keys(ligation_points)]), BenjaminiHochberg())
        for ((p, c), f) in zip(ligation_points, fdr)
            f < max_fdr || continue
            for pl in interact.bpstats[p][4]:interact.bpstats[p][5]
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

function aggregation_plot!(ax::Axis, interact::Interactions, n::String, max_fdr::Float64)
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
        if st == '+'
            allpositions .-= cds
        else
            allpositions = reverse(-1.0 .* (allpositions .- cds))
            reverse!(allcounts)
        end
        allpositions .+= 1
        band!(ax, allpositions, zeros(length(allcounts)), allcounts, label=label, color=color)
    end
end

function seed_correlation_plot!(ax::Axis, interact::Interactions, n::String, range1::UnitRange, range2::UnitRange, hn1::String, hn2::String, max_fdr::Float64)
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
        scatter!(ax, log.(collect(values(region1)) .+ 1), log.(collect(values(region2)) .+ 1), label=label, color=color)
        #density!(ax, log.(collect(values(region1)) .+ 1), color = color, bandwidth = 1.0)
        #density!(ax, log.(collect(values(region2)) .+ 1), color = color, bandwidth = 1.0, direction=:y)

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

function plot_figure_5(assets_folder::String, interact::Interactions)

    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 800*resfactor))

    Label(fig[1,1, TopLeft()], "A", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(fig[1,2, TopLeft()], "B", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    Label(fig[2,1, TopLeft()], "C", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,2, TopLeft()], "D", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    colors = reverse(ColorSchemes.tofino10.colors)[[4,3,2,1]]
    #colors = ColorSchemes.berlin25.colors[[7,5,3,1]]
    #colors = ColorSchemes.berlin25.colors[[3,5,7,9]]
    #colors = Makie.wong_colors()
    #colors = ColorSchemes.grayC10.colors[[2,4,6,8]]

    ax3 = Axis(fig[1, 1], title="FarS basepairing", ylabel="count", xlabel="position")
    aggregation_plot!(ax3, interact, "FarS", 0.25)
    band!(ax3, [30., 40.], [-200., -200.], [0., 0.], color=colors[3], label="seed 1")
    band!(ax3, [53., 60.], [-200., -200.], [0., 0.], color=colors[4], label="seed 2")
    ylims!(ax3, (-200, 4000))
    axislegend(ax3, position=:rt)


    ax_seed_corr = Axis(fig[1, 2], title="FarS seeds correlation", ylabel="log count in seed 2", xlabel="log count in seed 1")
    seed_correlation_plot!(ax_seed_corr, interact, "FarS", 30:40, 53:60, "VC1043", "VCA0848", 0.25)
    Legend(fig[1,3], ax_seed_corr)

    img4 = rotr90(load(joinpath(assets_folder, "basepairing.png")))
    ax4 = Axis(fig[2, 1], title="basepairing predictions",
        leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    hidedecorations!(ax4)
    image!(ax4, img4, aspect = DataAspect())

    df_plate = DataFrame(CSV.File(joinpath(assets_folder, "platereader2.csv")))

    ax5 = Axis(fig[2,2], title = "FarS reporter assay", ylabel="relative fluorescence [AU]", xticks = (1:2, [rich("vc1043"; font=:italic), rich("vca0848"; font=:italic)]))

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

    save( "figure_5.svg", fig)
    save( "figure_5.png", fig, px_per_unit = 2)

end
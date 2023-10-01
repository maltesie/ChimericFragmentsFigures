using CairoMakie, FileIO, MultipleTesting, JLD2, DataFrames, CSV, StatsBase, MultipleTesting, ColorSchemes, GeometryBasics, LaTeXStrings

struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int, Dict{Tuple{Int,Int},Int}, Dict{Tuple{Int,Int},Int}}}
    bpstats::Dict{Tuple{Int,Int}, Tuple{Float64, Int64, Int64, Int64, Int64, Float64}}
    multichimeras::Dict{Vector{Int}, Int}
    replicate_ids::Vector{Symbol}
    counts::Dict{Symbol,Vector{Int}}
end

Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

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

function plot_seed_regions()

    fig = Figure(resolution=(1200, 750))

    ga = fig[2, 1] = GridLayout()
    Label(ga[1,1, TopLeft()], "D", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,2, TopLeft()], "E", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,3, TopLeft()], "F", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    gb = fig[1, 1] = GridLayout()
    Label(gb[1,1, TopLeft()], "A", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(gb[1,3, TopLeft()], "B", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(gb[1,6, TopLeft()], "C", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    interact = Interactions("/home/abc/Workspace/ChimericFragmentsProjects/rilseq_vibrio_hfq_with_trimming/replicates_correlation/results_12_17_1-00/jld/hfq_hcd.jld2")

    #ax1 = Axis(ga[1, 1], title="OppZ", ylabel="count", xlabel="position")
    #aggregation_plot!(ax1, interact, "OppZ", 0.25)

    #ax2 = Axis(ga[1, 2], title="VC2770", ylabel="count", xlabel="position")
    #aggregation_plot!(ax2, interact, "VC2770", 0.25)

    ax3 = Axis(ga[1, 1], title="FarS basepairing", ylabel="count", xlabel="position")
    aggregation_plot!(ax3, interact, "FarS", 0.25)
    band!(ax3, [0., 0.], [0., 0.], [0., 0.], color=RGBAf(0.529, 0.721, 0.431, 1.0), label="seed 1")
    band!(ax3, [0., 0.], [0., 0.], [0., 0.], color=RGBAf(0.85, 0.89, 0.6, 1.0), label="seed 2")
    axislegend(ax3, position=:rt)
    #Legend(ga[1,4], ax3)

    img4 = rotr90(load(joinpath(@__DIR__, "basepairing.png")))
    ax4 = Axis(ga[1, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    hidedecorations!(ax4)
    image!(ax4, img4, aspect = DataAspect())

    colors = reverse(ColorSchemes.tofino10.colors)[[4,3,2,1]]#[[2,3,4,5]]#[[4,3,2,1]]
    df_plate = DataFrame(CSV.File(joinpath(@__DIR__, "platereader2.csv")))

    ax5 = Axis(ga[1,3], title = "FarS reporter assay", ylabel="relative fluorescence [AU]", xticks = (1:2, [L"\textit{vc1043}", L"\textit{vca0848}"]))

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
    Legend(ga[1,4], elements, labels)
    #axislegend(ax5, elements, labels, position=:rt)

    df_old = DataFrame(CSV.File(joinpath(@__DIR__, "hfq_lcd_old.csv")))
    df_new = DataFrame(CSV.File(joinpath(@__DIR__, "hfq_lcd_new.csv")))

    tested_pairs = [
        ("VCA0946", "Spot42"),
        ("aphA", "Spot42"),
        ("cyaA", "Spot42"),
        ("VC0972", "Spot42"),
        ("Spot42", "VC0972"),
        ("VCA0181", "Spot42"),
        ("VC1043", "Spot42"),
        ("VC0966", "Spot42"),
        ("VC2030", "Spot42"),
        ("VCA0166", "Spot42"),
        ("VC1904", "Spot42"),
        ("VC2602", "Spot42"),
    ]
    old_pairs = [(row["name1"], row["name2"]) for row in eachrow(df_old)]

    tested_index = collect(findall(row->(row["name1"], row["name2"]) in tested_pairs, eachrow(df_new)))
    old_index = collect(findall(row->(row["name1"], row["name2"]) in old_pairs, eachrow(df_new)))
    new_index = vcat(tested_index, [i for i in 1:nrow(df_new) if !((i in old_index) | (i in tested_index))])

    df_total_old = DataFrame(CSV.File(joinpath(@__DIR__, "table_s1_lcd.csv")))
    df_total_new = DataFrame(CSV.File(joinpath(@__DIR__, "interactions_hfq_lcd_new.csv")))

    types = (["sRNA"], ["IGR"], ["CDS", "5UTR", "3UTR", "CDS_UTRS"])

    counter = zeros(Float64, 12)
    for row in eachrow(df_total_old)
        (ismissing(row.type1) || ismissing(row.type2)) && continue
        for (i, t) in enumerate(types)
            counter[i] += Int(row.type1 in t)
            counter[i] += Int(row.type2 in t)
        end
    end
    for row in eachrow(df_total_new)
        for (i, t) in enumerate(types)
            if row.nb_ints >= 20
                counter[i+6] += Int(row.type1 in t)
                counter[i+6] += Int(row.type2 in t)
                if row.fisher_fdr <= 0.05
                    counter[i+3] += Int(row.type1 in t)
                    counter[i+3] += Int(row.type2 in t)
                end
            end
            #counter[i+9] += Int(row.type1 in t)
            #counter[i+9] += Int(row.type2 in t)
            if row.bp_fdr <= 0.25
                counter[i+9] += Int(row.type1 in t)
                counter[i+9] += Int(row.type2 in t)
            end
        end
    end
    counter ./= 2
    ax1 = Axis(gb[1,1:2], ylabel="count", title="annotation types", xticks = (1:4, ["previous", "improved", "no FDR, reads>=20", "FDR<0.25, reads>=3"]), xticklabelrotation = pi/8)
    colors = Makie.wong_colors()
    groups = [1,2,3,1,2,3,1,2,3,1,2,3]

    barplot!(ax1, vcat(fill(1, 3), fill(2, 3), fill(3, 3), fill(4, 3)), counter, stack=groups, color=colors[groups])
    labels = ["sRNA", "IGR", "mRNA"]
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    axislegend(ax1, elements, labels, position=:lt)

    axb2 = Axis(gb[1,3:5], ylabel="-log10(complementarity FDR)", xlabel="log10(reads count)", title="Spot 42 interactions")

    scatter!(axb2, log10.(df_new[new_index, :nb_ints]), -1 .* log10.(df_new[new_index, :bp_fdr]),
        label="Fisher FDR > 0.05", color=RGBAf(1.0, 0.5, 0., 0.7), markersize=9)

    scatter!(axb2, log10.(df_new[old_index, :nb_ints]), -1 .* log10.(df_new[old_index, :bp_fdr]),
        label="Fisher FDR <= 0.05", color=RGBAf(0.4, 0.25, 0.8, 0.7), markersize=9)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.8))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.4))
    scatter!(axb2, log10.(df_new[tested_index, :nb_ints]), -1 .* log10.(df_new[tested_index, :bp_fdr]), label="picked for validation", marker=Polygon(p_big, [p_small]), markersize=7, color=RGBAf(0.8, 0.1, 0.0, 1.0)) # color=RGBAf(0.2, 0.5, 0.4, 0.9)
    #Legend(fig[1,3], ax2)
    axislegend(axb2, position=:lt)

    df_plate = DataFrame(CSV.File(joinpath(@__DIR__, "platereader.csv")))
    xlabels_latex = [L"\textit{%$cl}" for cl in df_plate.name]
    axb4 = Axis(gb[1,6:8], title="Spot 42 reporter assay", ylabel="relative fluorescence [AU]", xticks = (1:nrow(df_plate), xlabels_latex), xticklabelrotation = pi/4)
    groups = vcat(fill(1,nrow(df_plate)), fill(2,nrow(df_plate)))
    ctrl = Matrix(df_plate[:, [:ctrl1, :ctrl2, :ctrl3]])
    mean_ctrl = vec(mean(ctrl; dims=2))
    ctrl ./= mean_ctrl
    sd_ctrl = vec(std(ctrl; dims=2))
    spot = Matrix(df_plate[:, [:exp1, :exp2, :exp3]])
    spot ./= mean_ctrl
    mean_spot = vec(mean(spot; dims=2))
    sd_spot = vec(std(spot; dims=2))
    mean_ctrl = ones(nrow(df_plate))
    barplot!(axb4, vcat(1:nrow(df_plate), 1:nrow(df_plate)), vcat(mean_ctrl, mean_spot), dodge=groups, color=colors[groups])
    scatter!(axb4, collect(1:nrow(df_plate)) .- 0.35, vec(ctrl[:, 1]), color="black", markersize=5)
    scatter!(axb4, collect(1:nrow(df_plate)) .- 0.2, vec(ctrl[:, 2]), color="black", markersize=5)
    scatter!(axb4, collect(1:nrow(df_plate)) .- 0.05, vec(ctrl[:, 3]), color="black", markersize=5)
    errorbars!(collect(1:nrow(df_plate)) .- 0.2, mean_ctrl, sd_ctrl, whiskerwidth=5)
    scatter!(axb4, collect(1:nrow(df_plate)) .+ 0.35, vec(spot[:, 1]), color="black", markersize=5)
    scatter!(axb4, collect(1:nrow(df_plate)) .+ 0.2, vec(spot[:, 2]), color="black", markersize=5)
    scatter!(axb4, collect(1:nrow(df_plate)) .+ 0.05, vec(spot[:, 3]), color="black", markersize=5)
    errorbars!(collect(1:nrow(df_plate)) .+ 0.2, mean_spot, sd_spot, whiskerwidth=5)
    labels = ["control", "Spot 42"]
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    #Legend(fig[2,3], elements, labels)
    axislegend(axb4, elements, labels, position=:lt)

    save(joinpath(@__DIR__, "seed_regions_and_new_vs_old.svg"), fig)
end
plot_seed_regions()
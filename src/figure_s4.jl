using RNASeqTools, DataFrames, CairoMakie, GeometryBasics, FileIO, JLD2, MultipleTesting, StatsBase, GenomicFeatures, CSV, HypothesisTests

struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int, Dict{Tuple{Int,Int},Int}, Dict{Tuple{Int,Int},Int}}}
    bpstats::Dict{Tuple{Int,Int}, Tuple{Float64, Int64, Int64, Int64, Int64, Float64, Int}}
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
function coverage_plot!(ax::Axis, coverage::Coverage, interval::Interval)
    vals = values(coverage, interval)
    strand(interval) == '-' && reverse!(vals)
    band!(ax, collect(0:(rightposition(interval)-leftposition(interval))), zeros(rightposition(interval)-leftposition(interval)+1), vals)
end

function new_regulators_plot()
    resfactor = 1.0
    fig = Figure(resolution=(1200*resfactor, 1200*resfactor))

    interact = Interactions("/home/abc/Workspace/ChimericFragmentsProjects/rilseq_vibrio_hfq_with_trimming/results/jld/hfq_lcd.jld2")

    ax_graph = Axis(fig[1,1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(@__DIR__, "hfq_lcd_igr_vcr069_vc1803.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())


    ga = fig[2:3,1] = GridLayout()
    ax1_1 = Axis(ga[1, 1], ylabel="count", title="basepairing predictions")

    limits_low, limits_high = 1530, 1630

    aggregation_plot!(ax1_1, interact, "Vcr069:VC1803", 0.9)
    axislegend(ax1_1)
    hidexdecorations!(ax1_1, grid = false)
    hidespines!(ax1_1, :l, :t, :r)
    xlims!(ax1_1, limits_low, limits_high)

    idx = findfirst(interact.nodes.name .== "Vcr069:VC1803")
    ref = interact.nodes.ref[idx]
    st = interact.nodes.strand[idx]
    interval = Interval(ref, interact.nodes.left[idx], interact.nodes.right[idx], st)
    coverage_tex = Coverage(joinpath(@__DIR__, "..", "new_regulators", "trimmed_tex_01_1_forward.bw"), joinpath(@__DIR__, "..", "new_regulators", "trimmed_tex_01_1_reverse.bw"))
    coverage_term = Coverage(joinpath(@__DIR__, "..", "new_regulators",  "trimmed_12898-mock-rep1_S33_R1_001_forward.bw"), joinpath(@__DIR__, "..", "new_regulators", "trimmed_12898-mock-rep1_S33_R1_001_reverse.bw"))

    ax1_2 = Axis(ga[2,1], ylabel="count", title="dRNA-seq")
    coverage_plot!(ax1_2, coverage_tex, interval)
    hidexdecorations!(ax1_2, grid = false)
    hidespines!(ax1_2, :l, :t, :r)
    xlims!(ax1_2, limits_low, limits_high)

    ax1_3 = Axis(ga[3,1], ylabel="count", xlabel="position in IGR", title="TERM-seq")
    coverage_plot!(ax1_3, coverage_term, interval)
    hidespines!(ax1_3, :l, :t, :r)
    xlims!(ax1_3, limits_low, limits_high)

    linkxaxes!(ax1_1, ax1_2, ax1_3)

    img4 = rotr90(load(joinpath(@__DIR__, "northern_2.png")))
    ax4 = Axis(fig[1:3, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    hidedecorations!(ax4)
    image!(ax4, img4, aspect = DataAspect())

    Label(fig[1,1, TopLeft()], "A", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,1, TopLeft()], "B", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[1,2, TopLeft()], "C", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)


    save(joinpath(@__DIR__, "figure_si_S4.svg"), fig)

end
new_regulators_plot()
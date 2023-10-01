using RNASeqTools, BioAlignments, BioSequences, CairoMakie, DataFrames, JLD2, MultipleTesting, FileIO, Combinatorics, StatsBase

randstrand() = rand() > 0.5 ? '+' : '-'

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

score_bp(paln::PairwiseAlignmentResult, shift_weight::Float64) = BioAlignments.score(paln) - (shift_weight * abs(paln.aln.a.aln.anchors[end].seqpos - paln.aln.a.aln.anchors[end].refpos))
function random_bp_alignments_plot(nseqs::Int, check_interaction_distances::Tuple{Int,Int}, bp_parameters::NTuple{6,Int}, params::Vector{Tuple{Int, Int}})
    genome = Genome(joinpath(@__DIR__, "..", "bwa_mem2", "mg1655.fna"))
    #genome.seq = randdnaseq(length(genome.seq))
    #genome = Genome(randdnaseq(4000000), Dict("test"=>1:4000000))
    complement_genome = Genome(complement(genome.seq), genome.chroms)
    reverse_genome = Genome(copy(genome.seq), genome.chroms)
    for (_, seq) in reverse_genome
        reverse!(seq)
    end
    reverse_complement_genome = Genome(complement(genome.seq), genome.chroms)
    for (_, seq) in reverse_complement_genome
        reverse!(seq)
    end

    scores = Dict((DNA_A, DNA_T)=>bp_parameters[1], (DNA_T, DNA_A)=>bp_parameters[1],
                    (DNA_C, DNA_G)=>bp_parameters[2], (DNA_G, DNA_C)=>bp_parameters[2],
                    (DNA_G, DNA_T)=>bp_parameters[3], (DNA_T, DNA_G)=>bp_parameters[3])
    model = AffineGapScoreModel(SubstitutionMatrix(scores; default_match=-1*bp_parameters[4], default_mismatch=-1*bp_parameters[4]);
        gap_open=-1*bp_parameters[5], gap_extend=-1*bp_parameters[6])

    #model = AffineGapScoreModel(match=1,mismatch=-1,gap_open=-1,gap_extend=-1)
    #model = CostModel(match=-1,mismatch=1,insertion=1,deletion=1)
    #model = CostModel(SubstitutionMatrix(scores; default_match=bp_parameters[4], default_mismatch=bp_parameters[4]), bp_parameters[5], bp_parameters[6])
    left1, left2, right1, right2, rands = zeros(Int, nseqs), zeros(Int, nseqs), zeros(Int, nseqs), zeros(Int, nseqs), zeros(Float64, nseqs)

    bins= 1:(check_interaction_distances[1] - check_interaction_distances[2])+1
    #distcut = 5
    #ccc = 0
    #ops = Dict()
    for (i, (i1::Int, i2::Int)) in enumerate(eachrow(rand(1:length(genome.seq), (nseqs, 2))))

        strand1, strand2 = randstrand(), randstrand()
        ref1, ref2 = rand(keys(genome.chroms)), rand(keys(genome.chroms))
        l1, l2 = length(genome.chroms[ref1]), length(genome.chroms[ref2])
        s1 = if strand1=='+'
            view(genome[ref1], clamp(i1-check_interaction_distances[1]+1, 1, l1):clamp(i1-check_interaction_distances[2], 1, l1))
        else
            c1, c2 = l1+1-clamp(i1+check_interaction_distances[2], 1, l1), l1+1-clamp(i1+check_interaction_distances[1]-1, 1, l1)
            view(reverse_complement_genome[ref1], c2:c1)
        end
        s2 = if strand2=='-'
            view(complement_genome[ref2], clamp(i2-check_interaction_distances[1]+1, 1, l2):clamp(i2-check_interaction_distances[2], 1, l2))
        else
            c1, c2 = l2+1-clamp(i2+check_interaction_distances[2], 1, l2), l2+1-clamp(i2+check_interaction_distances[1]-1, 1, l2)
            view(reverse_genome[ref2], c2:c1)
        end

        s1 = randdnaseq(last(bins)-1)
        s2 = randdnaseq(last(bins)-1)

        paln = pairalign(LocalAlignment(), s2, s1, model)
        #println(BioAlignments.print_pairwise_alignment(stdout,paln.aln))
        rands[i] = score_bp(paln, 1.0)
        left1[i] = paln.aln.a.aln.anchors[1].seqpos + 1
        right1[i] = paln.aln.a.aln.anchors[end].seqpos
        left2[i] = paln.aln.a.aln.anchors[1].refpos + 1
        right2[i] = paln.aln.a.aln.anchors[end].refpos
        #left1[i] = paln.aln.a.aln.anchors[findfirst(x -> ismatchop(x.op), paln.aln.a.aln.anchors)-1].seqpos+1
        #right1[i] = paln.aln.a.aln.anchors[findlast(x -> ismatchop(x.op), paln.aln.a.aln.anchors)].seqpos
        #left2[i] = paln.aln.a.aln.anchors[findfirst(x -> ismatchop(x.op), paln.aln.a.aln.anchors)-1].refpos+1
        #right2[i] = paln.aln.a.aln.anchors[findlast(x -> ismatchop(x.op), paln.aln.a.aln.anchors)].refpos

        #ccc += (abs(right1[i] - right2[i]) <= distcut)
        #paln.aln.a.aln.anchors[end].op in keys(ops) ? ops[paln.aln.a.aln.anchors[end].op] += 1 : ops[paln.aln.a.aln.anchors[end].op] = 1
    end
    #println(ccc / nseqs)

    fig = Figure(resolution=(1400, 900))

    ax1 = Axis(fig[2,1], title="random, left", xlabel="postion in alignment", ylabel="frequency")
    h11 = hist!(ax1, left1, label="RNA1, left", bins=bins, normalization=:probability)
    h12 = hist!(ax1, left2, label="RNA2, left", bins=bins, normalization=:probability)

    ax2 = Axis(fig[2,2], title="random, right", xlabel="postion in alignment", ylabel="frequency")
    h21 = hist!(ax2, right1, label="RNA1, right", bins=bins, color=RGBAf(0.9, 0.2, 0.2, 0.6), normalization=:probability)
    h22 = hist!(ax2, right2, label="RNA2, right", bins=bins, color=RGBAf(0.5, 0.3, 0.9, 0.6), normalization=:probability)

    linkyaxes!(ax1, ax2)
    Label(fig[2,1, TopLeft()], "D", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    interact = Interactions(joinpath(@__DIR__, "..", "..", "..", "ChimericFragmentsProjects", "rilseq_vibrio_hfq_with_trimming", "results", "jld","hfq_lcd.jld2"))

    resize!.((left1, left2, right1, right2), length(interact.bpstats))
    fdrv = adjust(PValues(first.(values(interact.bpstats))), BenjaminiHochberg())
    interacts = zeros(Float64, length(fdrv))
    fcut = 0.25
    #cccc = 0
    #norm = 0
    #rs = Int64[]
    #ls = Int64[]
    for (i, (pval, l1, r1, l2, r2, s)) in enumerate(values(interact.bpstats))
        left1[i] = l1
        left2[i] = l2
        right1[i] = r1
        right2[i] = r2
        interacts[i] = s
        #if pval > 0.25
        #    if abs(right1[i] - right2[i]) <= distcut
        #        cccc += 1
        #        push!(rs, r2)
        #        push!(ls, l2)
        #    end
        #    norm += 1
        #end
    end
    #println(cccc/norm)
    #println(norm)


    #randseq_model_ecdf = ecdf(rands)
    #interact_ecdf = ecdf(interacts)

    #max_score = max(maximum(interacts), maximum(rands))
    #randseq_model_pdf = diff(randseq_model_ecdf.(1:(max_score+1)))
    #interact_pdf = diff(interact_ecdf.(1:(max_score+1)))

    ax9 = Axis(fig[2,3], xlabel="postion in alignment", ylabel="frequency", title="total, left")
    h91 = hist!(ax9, left1, label="RNA1, left", bins=bins, normalization=:probability)
    h92 = hist!(ax9, left2, label="RNA2, left", bins=bins, normalization=:probability)

    ax10 = Axis(fig[2,4], xlabel="postion in alignment", ylabel="frequency", title="total, right")
    h101 = hist!(ax10, right1, label="RNA1, right", bins=bins, color=RGBAf(0.9, 0.2, 0.2, 0.6), normalization=:probability)
    h102 = hist!(ax10, right2, label="RNA2, right", bins=bins, color=RGBAf(0.5, 0.3, 0.9, 0.6), normalization=:probability)

    linkyaxes!(ax9, ax10)
    Label(fig[2,3, TopLeft()], "E", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    ax3 = Axis(fig[3,3], xlabel="postion in alignment", ylabel="frequency", title="significant, left")
    h31 = hist!(ax3, left1[fdrv .<= fcut], label="RNA1, left", bins=bins, normalization=:probability)
    h32 = hist!(ax3, left2[fdrv .<= fcut], label="RNA2, left", bins=bins, normalization=:probability)

    ax4 = Axis(fig[3,4], xlabel="postion in alignment", ylabel="frequency", title="significant, right")
    h41 = hist!(ax4, right1[fdrv .<= fcut], label="RNA1, right", bins=bins, color=RGBAf(0.9, 0.2, 0.2, 0.6), normalization=:probability)
    h42 = hist!(ax4, right2[fdrv .<= fcut], label="RNA2, right", bins=bins, color=RGBAf(0.5, 0.3, 0.9, 0.6), normalization=:probability)

    linkyaxes!(ax3, ax4)
    Label(fig[3,3, TopLeft()], "G", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    ax5 = Axis(fig[3,1], xlabel="postion in alignment", ylabel="frequency", title="unsignificant, left")
    h51 = hist!(ax5, left1[fdrv .> fcut], label="RNA1, left", bins=bins, normalization=:probability)
    h52 = hist!(ax5, left2[fdrv .> fcut], label="RNA2, left", bins=bins, normalization=:probability)

    ax6 = Axis(fig[3,2], xlabel="postion in alignment", ylabel="frequency", title="unsignificant, right")
    h61 = hist!(ax6, right1[fdrv .> fcut], label="RNA1, right", bins=bins, color=RGBAf(0.9, 0.2, 0.2, 0.6), normalization=:probability)
    h62 = hist!(ax6, right2[fdrv .> fcut], label="RNA2, right", bins=bins, color=RGBAf(0.5, 0.3, 0.9, 0.6), normalization=:probability)

    linkyaxes!(ax5, ax6)
    Label(fig[3,1, TopLeft()], "F", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    Legend(fig[2:3,5], [h31, h32, h41, h42], ["RNA1, left", "RNA2, left", "RNA1, right", "RNA2, right"])

    ga = fig[1, 1:5]= GridLayout()

    significant_score = minimum(interacts[fdrv .<= fcut])
    ax7 = Axis(ga[1, 2], title="density of complementarity scores", xlabel="score", ylabel="density")
    density!(ax7, rands, color=(:red, 0.3), label="random model")
    density!(ax7, interacts, color=(:green, 0.3), label="ligation points")
    vlines!(ax7, [significant_score], color = :blue, label="FDR = $(fcut)")
    axislegend(ax7)
    Label(ga[1,2, TopLeft()], "B", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    #ax9 = Axis(ga[1, 3], title="density of alignment scores", xlabel="score", ylabel="density")
    #density!(ax9, rands, color=(:red, 0.3), label="random model")
    #density!(ax9, interacts, color=(:green, 0.3), label="ligation points")
    #vlines!(ax9, [significant_score], color = :blue, label="FDR = $(fcut)")
    #axislegend(ax9)


    interact_path = "/home/abc/Workspace/ChimericFragmentsProjects/rilseq_vibrio_hfq_with_trimming/replicates_correlation/"
    interact_file = "hfq_lcd.jld2"
    colors = ("Brown", "Coral", "BlueViolet", "DarkGreen")
    pcuts = [0.01, 0.05, 0.1, 0.25, 0.5, 1.0]
    ax_cor = Axis(ga[1, 3], ylabel="rank correlation", xlabel="complementarity FDR cutoff", title="RIL-seq replicate correlation",
        xticks=(1:length(pcuts)+1, [["$(round(pc, digits=2))" for pc in pcuts]..., "all"]))

    for ((se, ms), color) in zip(params, colors)
        folder = "results_$(se)_$(ms)_1-00"
        l = "$(se) | $(ms)"
        fp = joinpath(interact_path, folder, "jld", interact_file)
        interact = Interactions(fp)
        filter!(:nb_ints => x -> x>3, interact.edges)
        corr_mean = zeros(length(pcuts)+1)
        corr_sd = zeros(length(pcuts)+1)
        counts = zeros(length(pcuts)+1)
        subcounts = zeros(length(pcuts))
        corr_top_mean = zeros(length(pcuts))

        for (i, pcut) in enumerate(pcuts)
            pindex = interact.edges.bp_fdr .<= pcut
            corr = [corspearman(interact.edges[pindex, p1], interact.edges[pindex, p2]) for (p1, p2) in combinations(interact.replicate_ids, 2)]
            subpindex = sort(sample(1:findlast(pindex), sum(pindex), replace=false))
            count_ints = Int(floor(nrow(interact.edges)*pcut))
            subcounts[i] = count_ints
            subpindex = 1:count_ints

            corr_top = [corspearman(interact.edges[subpindex, p1], interact.edges[subpindex, p2]) for (p1, p2) in combinations(interact.replicate_ids, 2)]
            corr_mean[i] = mean(corr)
            corr_sd[i] = std(corr)
            counts[i] = sum(pindex)
            corr_top_mean[i] = mean(corr_top)
        end

        corr = [corspearman(interact.edges[!, p1], interact.edges[!, p2]) for (p1, p2) in combinations(interact.replicate_ids, 2)]
        corr_mean[length(pcuts)+1] = mean(corr)
        corr_sd[length(pcuts)+1] = std(corr)
        counts[length(pcuts)+1] = nrow(interact.edges)
        scatter!(ax_cor, 1:(length(pcuts)+1), corr_mean, label=l, color=color)
    end
    Legend(ga[1,4], ax_cor, "seed | score")
    Label(ga[1,3, TopLeft()], "C", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    img = rotr90(load(joinpath(@__DIR__, "combined.png")))
    ax8 = Axis(ga[1, 1],
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false, aspect = DataAspect()
    )
    hidedecorations!(ax8)
    image!(ax8, img, aspect = DataAspect())
    Label(ga[1,1, TopLeft()], "A", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    save(joinpath(@__DIR__, "figure_bp_predictions.svg"), fig)
end

random_bp_alignments_plot(100000, (30,0), (4,5,0,7,8,3), [(12, 14), (15, 14), (12, 17),  (15, 18)])
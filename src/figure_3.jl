randstrand() = rand() > 0.5 ? '+' : '-'
score_bp(paln::PairwiseAlignmentResult, shift_weight::Float64) = BioAlignments.score(paln) - (shift_weight * abs(paln.aln.a.aln.anchors[end].seqpos - paln.aln.a.aln.anchors[end].refpos))
function plot_figure_3(assets_folder::String, interact::InteractionsNew, nseqs::Int, check_interaction_distances::Tuple{Int,Int}, bp_parameters::NTuple{6,Int}, params::Vector{Tuple{Int, Int}})
    genome = Genome(joinpath(assets_folder, "mg1655.fna"))
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

    left1, left2, right1, right2, rands = zeros(Int, nseqs), zeros(Int, nseqs), zeros(Int, nseqs), zeros(Int, nseqs), zeros(Float64, nseqs)

    bins= 1:(check_interaction_distances[1] - check_interaction_distances[2])+1

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
        rands[i] = score_bp(paln, 1.0)
        left1[i] = paln.aln.a.aln.anchors[1].seqpos + 1
        right1[i] = paln.aln.a.aln.anchors[end].seqpos
        left2[i] = paln.aln.a.aln.anchors[1].refpos + 1
        right2[i] = paln.aln.a.aln.anchors[end].refpos
    end

    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 850*resfactor))

    ax1 = Axis(fig[2,1], title="random, left", xlabel="postion in alignment", ylabel="frequency")
    h12 = hist!(ax1, left2, label="RNA2, left", bins=bins, normalization=:probability, color=RGBAf(0.937, 0.333, 0.231, 0.6))
    h11 = hist!(ax1, left1, label="RNA1, left", bins=bins, normalization=:probability, color=RGBAf(0.388, 0.431, 0.98, 0.6))

    ax2 = Axis(fig[2,2], title="random, right", xlabel="postion in alignment", ylabel="frequency")
    h22 = hist!(ax2, right2, label="RNA2, right", bins=bins, normalization=:probability, color=RGBAf(0.937, 0.333, 0.231, 0.6))
    h21 = hist!(ax2, right1, label="RNA1, right", bins=bins, normalization=:probability, color=RGBAf(0.388, 0.431, 0.98, 0.6))

    #linkyaxes!(ax1, ax2)


    resize!.((left1, left2, right1, right2), length(interact.bpstats))
    fdrv = adjust(PValues(first.(values(interact.bpstats))), BenjaminiHochberg())
    interacts = zeros(Float64, length(fdrv))
    fcut = 0.25

    for (i, (pval, l1, r1, l2, r2, s)) in enumerate(values(interact.bpstats))
        left1[i] = l1
        left2[i] = l2
        right1[i] = r1
        right2[i] = r2
        interacts[i] = s
    end

    ax9 = Axis(fig[2,3], xlabel="postion in alignment", ylabel="frequency", title="total, left")
    h92 = hist!(ax9, left2, label="RNA2, left", bins=bins, normalization=:probability, color=RGBAf(0.937, 0.333, 0.231, 0.6))
    h91 = hist!(ax9, left1, label="RNA1, left", bins=bins, normalization=:probability, color=RGBAf(0.388, 0.431, 0.98, 0.6))

    ax10 = Axis(fig[2,4], xlabel="postion in alignment", ylabel="frequency", title="total, right")
    h102 = hist!(ax10, right2, label="RNA2, right", bins=bins, normalization=:probability, color=RGBAf(0.937, 0.333, 0.231, 0.6))
    h101 = hist!(ax10, right1, label="RNA1, right", bins=bins, normalization=:probability, color=RGBAf(0.388, 0.431, 0.98, 0.6))

    #linkyaxes!(ax9, ax10)


    ax3 = Axis(fig[3,3], xlabel="postion in alignment", ylabel="frequency", title="significant, left")
    h32 = hist!(ax3, left2[fdrv .<= fcut], label="RNA2, left", bins=bins, normalization=:probability, color=RGBAf(0.937, 0.333, 0.231, 0.6))
    h31 = hist!(ax3, left1[fdrv .<= fcut], label="RNA1, left", bins=bins, normalization=:probability, color=RGBAf(0.388, 0.431, 0.98, 0.6))

    ax4 = Axis(fig[3,4], xlabel="postion in alignment", ylabel="frequency", title="significant, right")
    h42 = hist!(ax4, right2[fdrv .<= fcut], label="RNA2, right", bins=bins, normalization=:probability, color=RGBAf(0.937, 0.333, 0.231, 0.6))
    h41 = hist!(ax4, right1[fdrv .<= fcut], label="RNA1, right", bins=bins, normalization=:probability, color=RGBAf(0.388, 0.431, 0.98, 0.6))

    #linkyaxes!(ax3, ax4)


    ax5 = Axis(fig[3,1], xlabel="postion in alignment", ylabel="frequency", title="unsignificant, left")
    h52 = hist!(ax5, left2[fdrv .> fcut], label="RNA2, left", bins=bins, normalization=:probability, color=RGBAf(0.937, 0.333, 0.231, 0.6))
    h51 = hist!(ax5, left1[fdrv .> fcut], label="RNA1, left", bins=bins, normalization=:probability, color=RGBAf(0.388, 0.431, 0.98, 0.6))

    ax6 = Axis(fig[3,2], xlabel="postion in alignment", ylabel="frequency", title="unsignificant, right")
    h62 = hist!(ax6, right2[fdrv .> fcut], label="RNA2, right", bins=bins, normalization=:probability, color=RGBAf(0.937, 0.333, 0.231, 0.6))
    h61 = hist!(ax6, right1[fdrv .> fcut], label="RNA1, right", bins=bins, normalization=:probability, color=RGBAf(0.388, 0.431, 0.98, 0.6))

    #linkyaxes!(ax5, ax6)

    linkyaxes!(ax1, ax2, ax3, ax4, ax5, ax6, ax9, ax10)

    Legend(fig[2:3,5], [h31, h32], ["RNA1", "RNA2"])

    ga = fig[1, 1:5]= GridLayout()

    significant_score = minimum(interacts[fdrv .<= fcut])
    ax7 = Axis(ga[1, 2], title="complementarity scores", xlabel="score", ylabel="density")
    density!(ax7, rands, color=(:red, 0.3), label="random model")
    density!(ax7, interacts, color=(:green, 0.3), label="ligation points")
    vlines!(ax7, [significant_score], color = :blue, label="FDR = $(fcut)")
    axislegend(ax7)


    colors = ("Brown", "Coral", "BlueViolet", "DarkGreen")
    pcuts = [0.01, 0.05, 0.1, 0.25, 0.5, 1.0]
    ax_cor = Axis(ga[1, 3], ylabel="rank correlation", xlabel="complementarity FDR cutoff", title="RIL-seq replicate correlation",
        xticks=(1:length(pcuts)+1, [["$(round(pc, digits=2))" for pc in pcuts]..., "all"]))
    replicate_ids = ["hfq_lcd_1", "hfq_lcd_2"]

    for ((se, ms), color) in zip(params, colors)
        l = "$(se) | $(ms)"
        fp = joinpath(assets_folder, "..", "figure_2", "csv_correlation", "hfq_lcd_$(se)_$(ms).csv")
        df = DataFrame(CSV.File(fp))
        corr_mean = zeros(length(pcuts)+1)
        corr_sd = zeros(length(pcuts)+1)
        counts = zeros(length(pcuts)+1)
        subcounts = zeros(length(pcuts))
        corr_top_mean = zeros(length(pcuts))

        for (i, pcut) in enumerate(pcuts)
            pindex = df.bp_fdr .<= pcut
            corr = [corspearman(df[pindex, p1], df[pindex, p2]) for (p1, p2) in combinations(replicate_ids, 2)]
            subpindex = sort(sample(1:findlast(pindex), sum(pindex), replace=false))
            count_ints = Int(floor(nrow(df)*pcut))
            subcounts[i] = count_ints
            subpindex = 1:count_ints

            corr_top = [corspearman(df[subpindex, p1], df[subpindex, p2]) for (p1, p2) in combinations(replicate_ids, 2)]
            corr_mean[i] = mean(corr)
            corr_sd[i] = std(corr)
            counts[i] = sum(pindex)
            corr_top_mean[i] = mean(corr_top)
        end

        corr = [corspearman(df[!, p1], df[!, p2]) for (p1, p2) in combinations(replicate_ids, 2)]
        corr_mean[length(pcuts)+1] = mean(corr)
        corr_sd[length(pcuts)+1] = std(corr)
        counts[length(pcuts)+1] = nrow(df)
        scatter!(ax_cor, 1:(length(pcuts)+1), corr_mean, label=l, color=color)
    end

    Legend(ga[1,4], ax_cor, "seed | score")

    img = rotr90(load(joinpath(assets_folder, "combined.png")))
    ax8 = Axis(ga[1, 1], title="basepairing scheme")
    hidedecorations!(ax8)
    image!(ax8, img, aspect = DataAspect())
    Label(ga[1,1, TopLeft()], "a", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,2, TopLeft()], "b", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,3, TopLeft()], "c", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    Label(fig[2,1, TopLeft()], "d", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,3, TopLeft()], "e", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,1, TopLeft()], "f", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,3, TopLeft()], "g", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    #titlelayout = GridLayout(fig[0, 1], halign = :left, tellwidth = false)
    #Label(titlelayout[1, 1], "Fig. 3", halign = :left, fontsize=30)

    save("figure_3.pdf", fig)
    save("figure_3.png", fig, px_per_unit = 2)
end
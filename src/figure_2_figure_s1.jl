function make_random_single_testseqs(genome::Genome, unique_set::Set{Int}; nseqs=1000000, len=30, nerr=0)
    newseq = randdnaseq(nseqs*len)
    seqnames = Vector{UInt}(undef, nseqs)
    ranges = Vector{UnitRange{Int}}(undef, nseqs)
    truep = rand(unique_set, nseqs)
    for (i, p) in enumerate(truep)
        p in unique_set || continue
        newp = (((i-1) * len)+1)
        newr = newp:(newp+len-1)
        news = view(newseq, newr)
        copyto!(news, view(genome.seq, p:(p+len-1)))
        for i in 1:nerr
            news[rand(1:len)] = rand((DNA_A, DNA_G, DNA_C, DNA_T))
        end
        ranges[i] = newr
        seqnames[i] = UInt(i)
    end
    return Sequences(newseq, seqnames, ranges), truep
end

function linear_unique(n::Int, unique_set::Set{Int}; repeat_each=100)
    is = zeros(Int, n*repeat_each)
    c = 1
    for i in 1:n
        while !(c in unique_set)
            c += 1
        end
        is[((i-1)*repeat_each+1):(i*repeat_each)] .= c
        c += 1
    end
    is
end

function make_linear_single_testseqs(genome::Genome, unique_set::Set{Int}; nseqs=1000000, len=30, nerr=0)
    newseq = randdnaseq(nseqs*len)
    seqnames = Vector{UInt}(undef, nseqs)
    ranges = Vector{UnitRange{Int}}(undef, nseqs)
    truep = linear_unique(nseqs, unique_set)
    for (i, p) in enumerate(truep)
        p in unique_set || continue
        newp = (((i-1) * len)+1)
        newr = newp:(newp+len-1)
        news = view(newseq, newr)
        copyto!(news, view(genome.seq, p:(p+len-1)))
        for i in 1:nerr
            news[rand(1:len)] = rand((DNA_A, DNA_G, DNA_C, DNA_T))
        end
        ranges[i] = newr
        seqnames[i] = UInt(i)
    end
    return Sequences(newseq, seqnames, ranges), truep
end

function make_random_linear_single_testseqs(genome::Genome, unique_set::Set{Int}; nseqs=1000000, minlen=15, maxlen=40, maxrandlen=5, maxerr=1, repeateach=100, onlyunique=true)
    newseq = randdnaseq(nseqs*(maxlen+maxrandlen))
    seqnames = Vector{UInt}(undef, nseqs)
    ranges = Vector{UnitRange{Int}}(undef, nseqs)
    truep = onlyunique ? linear_unique(Int(floor(nseqs/repeateach)), unique_set; repeat_each=repeateach) : [i for i in Int(floor(nseqs/repeateach)) for _ in 1:repeateach]
    newp = 1
    for (i, p) in enumerate(truep)
        p in unique_set || continue
        len=rand(minlen:maxlen)
        newr = newp:(newp+len-1)
        news = view(newseq, newr)
        copyto!(news, view(genome.seq, p:(p+len-1)))
        nerr=rand(0:maxerr)
        for i in 1:nerr
            news[rand(1:len-1)] = rand((DNA_A, DNA_G, DNA_C, DNA_T))
        end
        ranges[i] = newp:(newp+len-1+rand(0:maxrandlen))
        newp = last(ranges[i])+1
        seqnames[i] = UInt(i)
    end
    return Sequences(newseq, seqnames, ranges), truep
end

function make_chimeric_testseqs(genome::Genome, unique_set::Set{Int}; nseqs=1000000, len1=11, len2=19, nerr=0)
    newseq = randdnaseq(nseqs*(len1+len2))
    seqnames = Vector{UInt}(undef, nseqs)
    ranges = Vector{UnitRange{Int}}(undef, nseqs)
    truep = rand(unique_set, (nseqs, 2))
    for (i, (p1, p2)) in enumerate(eachrow(truep))
        newp = (((i-1) * (len1+len2))+1)
        newr1 = newp:(newp+len1-1)
        news1 = view(newseq, newr1)
        copyto!(news1, view(genome.seq, p1:(p1+len1-1)))
        newr2 = (newp+len1):(newp+len1+len2-1)
        news2 = view(newseq, newr2)
        copyto!(news2, view(genome.seq, p2:(p2+len2-1)))
        for i in 1:nerr
            news1[rand(1:len1-1)] = rand((DNA_A, DNA_G, DNA_C, DNA_T))
            news2[rand(1:len2-1)] = rand((DNA_A, DNA_G, DNA_C, DNA_T))
        end
        ranges[i] = newp:(newp+len1+len2-1)
        seqnames[i] = UInt(i)
    end
    return Sequences(newseq, seqnames, ranges), truep
end

function make_random_chimeric_testseqs(genome::Genome, unique_set::Set{Int}; nseqs=1000000, minlen=13, maxlen=40, maxrandlen=0, maxerr=1)
    newseq = randdnaseq(nseqs*2*maxlen)
    seqnames = Vector{UInt}(undef, nseqs)
    ranges = Vector{UnitRange{Int}}(undef, nseqs)
    truep = rand(unique_set, (nseqs, 2))
    newp = 1
    for (i, (p1, p2)) in enumerate(eachrow(truep))
        len1=rand(minlen:maxlen)
        len2=rand(minlen:maxlen)
        newr1 = newp:(newp+len1-1)
        news1 = view(newseq, newr1)
        copyto!(news1, view(genome.seq, p1:(p1+len1-1)))
        newr2 = (newp+len1):(newp+len1+len2-1)
        news2 = view(newseq, newr2)
        copyto!(news2, view(genome.seq, p2:(p2+len2-1)))
        nerr=rand(0:maxerr)
        for i in 1:nerr
            news1[rand(1:len1-1)] = rand((DNA_A, DNA_G, DNA_C, DNA_T))
            news2[rand(1:len2-1)] = rand((DNA_A, DNA_G, DNA_C, DNA_T))
        end
        ranges[i] = newp:(newp+len1+len2-1+rand(0:maxrandlen))
        newp = last(ranges[i])+1
        seqnames[i] = UInt(i)
    end
    return Sequences(newseq, seqnames, ranges), truep
end

function bwa_mem_tpr_fpr(s::Sequences, truep::Matrix{Int}, g::Genome, min_seedlens::Vector{Int}, min_scores::Vector{Int})
    nseqs = length(s)
    tptrans = Dict(hash("$i")=>i for i in 1:nseqs)
    fasta1_path = joinpath(@__DIR__, "test_chimeric_1.fasta.gz")
    write(fasta1_path, s)
    reverse_complement!(s)
    fasta2_path = joinpath(@__DIR__, "test_chimeric_2.fasta.gz")
    write(fasta2_path, s)
    tpr_matrix = zeros((length(min_seedlens), length(min_scores)))
    fpr_matrix = zeros((length(min_seedlens), length(min_scores)))
    foundids = Dict{Tuple{Int,Int}, Set{UInt}}()
    for (i,min_seedlen) in enumerate(min_seedlens), (j,min_score) in enumerate(min_scores)
        align_mem(PairedSingleTypeFiles([fasta1_path], [fasta2_path]), g; min_seed_len=min_seedlen, min_score=min_score, unpair_penalty=15,
            output_all_alignments=false, unpair_rescue=true, overwrite_existing=true)
        alns = AlignedReads(joinpath(@__DIR__, "test_chimeric.bam"); include_alternative_alignments=false)
        cc = zeros(Int, 2)
        found = Set{UInt}()
        for aln in alns
            i1, i2 = truep[tptrans[readid(aln)], :]
            f1, f2 = any(abs(leftposition(alnp) - i1) <= 8 for alnp in aln), any(abs(leftposition(alnp) - i2) <= 8 for alnp in aln)
            if ischimeric(aln; check_annotation=false)
                cc[1] += 1
                if (f1 & f2)
                    cc[2] += 1
                    push!(found, readid(aln))
                end
            end
        end
        tpr_matrix[i, j] = cc[2]/nseqs
        fpr_matrix[i, j] = (cc[1]-cc[2])/nseqs
        foundids[(i, j)] = found
    end
    unique_rates = similar(tpr_matrix)
    for (i,j) in keys(foundids)
        union_rest = foldl(union!, s for (k,s) in foundids if ((k[1] == i) && (k[2] != j)))
        unique_ids = setdiff(foundids[(i,j)], union_rest)
        unique_rates[i,j] = length(unique_ids)/nseqs
    end
    return tpr_matrix, fpr_matrix, unique_rates
end

function bwa_mem_tpr_fpr(fname_tpr::String, fname_fpr::String, fname_unique::String, s::Sequences, truep::Matrix{Int}, g::Genome, seedlens::Vector{Int}, scores::Vector{Int})
    if isfile(fname_tpr) && isfile(fname_fpr) && isfile(fname_unique)
        readdlm(fname_tpr, '\t', Float64, '\n'), readdlm(fname_fpr, '\t', Float64, '\n'), readdlm(fname_unique, '\t', Float64, '\n')
    else
        tpr_tmp, fpr_tmp, unique_tmp = bwa_mem_tpr_fpr(s, truep, g, seedlens, scores)
        open(fname_tpr, "w") do io
            writedlm(io, tpr_tmp)
        end
        open(fname_fpr, "w") do io
            writedlm(io, fpr_tmp)
        end
        open(fname_unique, "w") do io
            writedlm(io, unique_tmp)
        end
        tpr_tmp, fpr_tmp, unique_tmp
    end
end

function conditional_distance_distribution(s::Sequences, truep::Vector{Int}, genome_path::String, params::Vector{Tuple{Int,Int}})
    nseqs = length(s)
    tptrans = Dict(hash("$i")=>i for i in 1:nseqs)
    fasta1_path = joinpath(@__DIR__, "test_linear.fasta.gz")
    !isfile(fasta1_path) && write(fasta1_path, s)
    mm_pairss = []
    lens = Int[]
    for (min_seed, min_score) in params
        mm_pairs = Tuple{Int,Int}[]
        bam_path = joinpath(@__DIR__, "test_linear_$min_seed-$min_score.bam")
        !isfile(bam_path) && align_mem(fasta1_path, genome_path, bam_path; min_seed_len=min_seed, min_score=min_score, output_all_alignments=false)
        alns = AlignedReads(bam_path; include_secondary_alignments=true)
        push!(lens, length(alns))
        for aln in alns
            i1 = truep[tptrans[readid(aln)]]
            length(aln) == 1 || continue
            lp = leftposition(first(aln))
            d = abs(lp - i1)
            d > 8 && push!(mm_pairs, (i1,lp))
        end
        push!(mm_pairss, mm_pairs)
    end
    println("found $(length.(mm_pairss)) misalignments in $lens alignments")
    return mm_pairss, lens
end

function probabilities_distances(genome_path::String, unique_set::Set{Int}, params::Vector{Tuple{Int,Int}};
        nseqs=1000000, minlen=13, maxlen=40, maxrandlen=0, maxerr=1, maxd=10000, repeateach=20, nsamples=10000)
    g = Genome(genome_path)
    s, truep = make_random_linear_single_testseqs(g, unique_set; nseqs=nseqs, minlen=minlen, maxlen=maxlen, maxrandlen=maxrandlen, maxerr=maxerr, repeateach=repeateach)
    mm_pairss, lens = conditional_distance_distribution(s, truep, genome_path, params)

    rates = [
        begin
            count_similar = zeros(Float64, length(100:100:maxd))
            correct = first.(mm_pairs)
            mm = last.(mm_pairs)
            sortindex = sortperm(correct)
            correct = correct[sortindex]
            mm = mm[sortindex]
            for (ii,d) in enumerate(100:100:maxd)
                low = 1
                high = 1
                for i in eachindex(correct)
                    low = findnext(x -> correct[i]-x <= d, correct, low)
                    high = findnext(x -> x-correct[i] > d, correct, high)
                    high = isnothing(high) ? length(correct) : high - 1
                    (i in low:high) || throw(AssertionError("weird"))
                    close_mismatches = sum(abs.(view(mm, low:high) .- mm[i]) .<= d) - 1
                    count_similar[ii] += close_mismatches
                end
            end
            count_similar ./ len
        end for (mm_pairs, len) in zip(mm_pairss, lens)
    ]
    return rates
end

function unique_set(genome::Genome; k=12, dist_to_max=500)
    us = Set{Int}()
    hs = Dict{UInt,Int}()
    for i in 1:length(genome.seq)-dist_to_max
        h = hash(view(genome.seq, i:(i+k-1)))
        if h in keys(hs)
            hs[h] in us && delete!(us, hs[h])
        else
            push!(us, i)
            push!(hs, h=>i)
        end
    end
    return us
end

function plot_figure_2(assets_path::String, lens::Vector{Int}, lmin::Int, lmax::Int, lmaxadapter::Int, nerrmax::Int, nseqs::Int,
        seedlens::Vector{Int}, scores::Vector{Int}, maxd::Int, params::Vector{Pair{Tuple{Int,Int}, String}})
    genome_path = joinpath(assets_path, "rand_genome.fna")
    if !isfile(genome_path)
        genomeseq = randdnaseq(5000000)
        g = Genome(genomeseq, Dict("test"=>1:5000000))
        write(genome_path, g)
    end
    g = Genome(genome_path)

    fig = Figure(resolution=(1300, 750))
    ga = fig[1, 1] = GridLayout()
    Label(ga[1,1, TopLeft()], "A", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(ga[1,3, TopLeft()], "B", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(ga[1,4, TopLeft()], "C", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    gb = fig[2, 1] = GridLayout()
    Label(gb[1,1, TopLeft()], "D", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(gb[1,2, TopLeft()], "E", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(gb[1,3, TopLeft()], "F", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)

    ax = Axis(ga[1,1], xlabel="FPR", ylabel="TPR", title="synthetic benchmarks")
    mkpath(joinpath(assets_path, "csv"))
    us = unique_set(g)

    fig_si = Figure(resolution=(1100, 1000))
    gc = fig_si[3, 1:4] = GridLayout()
    Label(fig_si[1,1, TopLeft()], "A", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(fig_si[1,2, TopLeft()], "B", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(fig_si[1,3, TopLeft()], "C", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(fig_si[2,1, TopLeft()], "D", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(fig_si[2,2, TopLeft()], "E", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(fig_si[2,3, TopLeft()], "F", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(gc[1,1, TopLeft()], "G", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(gc[1,2, TopLeft()], "H", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    Label(gc[1,3, TopLeft()], "I", fontsize = 26,font = :bold,padding = (0, 10, 10, 0), halign = :right)
    n_plots = 1
    for (j,nerr) in enumerate(0:nerrmax), (i,l) in enumerate(lens)
        s, truep = make_chimeric_testseqs(g, us; nseqs=nseqs, len1=l, len2=l, nerr=nerr)
        name = "length=$l, errors=$nerr"
        label = "$l | $nerr"
        fname_tpr = joinpath(assets_path, "csv", name * "_tpr.csv")
        fname_fpr = joinpath(assets_path, "csv", name * "_fpr.csv")
        fname_found = joinpath(assets_path, "csv", name * "_unique.csv")
        tpr, fpr, _ = bwa_mem_tpr_fpr(fname_tpr, fname_fpr, fname_found, s, truep, g, seedlens, scores)
        sorted_points = sort(collect(zip(vec(fpr), vec(tpr))), by=x->x[2])
        scatter!(ax, [first(p) for p in sorted_points], [last(p) for p in sorted_points], label=label, mode="markers", alpha=0.6)

        if nerr > 0
            clims = (0.0, 1.0)
            ax_si_tpr = Axis(fig_si[1, n_plots], yticks=(1:length(scores), string.(scores)), xticks = (1:length(seedlens), string.(seedlens)),
                ylabel="minimum alignment score", xlabel="seed length", title="TPR, length=$l")
            heatmap!(ax_si_tpr, tpr, colorrange=clims)
            for i in 1:length(seedlens), j in 1:length(scores)
                text!(ax_si_tpr, string(round(tpr[i,j], digits=2)), position = (i,j), align = (:center, :center), color = tpr[i,j] < 0.5 ? :white : :black)
            end
            ax_si_fpr = Axis(fig_si[2, n_plots], yticks=(1:length(scores), string.(scores)), xticks = (1:length(seedlens), string.(seedlens)),
                ylabel="minimum alignment score", xlabel="seed length", title="FPR, length=$l")
            heatmap!(ax_si_fpr, fpr, colorrange=clims)
            for i in 1:length(seedlens), j in 1:length(scores)
                text!(ax_si_fpr, string(round(fpr[i,j], digits=3)), position = (i,j), align = (:center, :center), color = fpr[i,j] < 0.5 ? :white : :black)
            end

            n_plots += 1
        end
    end

    Colorbar(fig_si[1:2,4], limits=(0.0, 1.0))

    s, truep = make_random_chimeric_testseqs(g, us; minlen=lmin, maxlen=lmax, maxrandlen=lmaxadapter, nseqs=nseqs)
    name = "length=$lmin-$lmax, errors=0-$nerrmax"
    label = "$lmin-$lmax | 0-$nerrmax"
    fname_tpr = joinpath(assets_path, "csv", name * "_tpr.csv")
    fname_fpr = joinpath(assets_path, "csv", name * "_fpr.csv")
    fname_found = joinpath(assets_path, "csv", name * "_unique.csv")
    tpr, fpr, _ = bwa_mem_tpr_fpr(fname_tpr, fname_fpr, fname_found, s, truep, g, seedlens, scores)
    clims = (min(minimum(tpr), minimum(fpr)), max(maximum(tpr), maximum(fpr)))
    noisex = (rand()-0.5)*0.0003
    noisey = (rand()-0.5)*0.0003
    sorted_points = sort(collect(zip(vec(fpr) .+ noisex, vec(tpr) .+ noisey)), by=x->x[2])
    scatter!(ax, [first(p) for p in sorted_points], [last(p) for p in sorted_points], label=label, mode="markers", alpha=0.6)
    Legend(ga[1,2], ax, "length | errors")

    interact_path = "/home/abc/Workspace/ChimericFragmentsProjects/rilseq_vibrio_hfq_with_trimming/replicates_correlation/"
    interact_file = "hfq_lcd.jld2"
    colors = ("Brown", "Coral", "BlueViolet", "DarkGreen")
    pcuts = [0.01, 0.05, 0.1, 0.25, 0.5, 1.0]
    ax_cor = Axis(gc[1,3], ylabel="Pearson correlation", xlabel="complementarity FDR cutoff", title="LCD replicate correlation", xticks=(1:length(pcuts)+1, [["$(round(pc, digits=2))" for pc in pcuts]..., "all"]))
    ax_count = Axis(gc[1,2], ylabel="median of read counts", xlabel="complementarity FDR cutoff", title="LCD reads per interaction", xticks=(1:length(pcuts)+1, [["$(round(pc, digits=2))" for pc in pcuts]..., "all"]), yscale=log10)
    ax_top = Axis(gb[1,3], ylabel="rank correlation", xlabel="top fraction of dataset", title="LCD replicate correlation", xticks=(1:length(pcuts), ["$(round(pc, digits=2))" for pc in pcuts]))
    ax_ints_count = Axis(gc[1,1], ylabel="median of read counts", xlabel="top fraction of dataset", title="LCD reads per interaction", xticks=(1:length(pcuts), ["$(round(pc, digits=2))" for pc in pcuts]), yscale=log10)
    max_count = 0

    for (((se, ms), label), color) in zip(params, colors)
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
        sorted_ints = sort(interact.edges.nb_ints; rev=true)
        for (i, pcut) in enumerate(pcuts)
            pindex = interact.edges.bp_fdr .<= pcut
            corr = [cor(interact.edges[pindex, p1], interact.edges[pindex, p2]) for (p1, p2) in combinations(interact.replicate_ids, 2)]
            subpindex = sort(sample(1:findlast(pindex), sum(pindex), replace=false))
            count_ints = Int(floor(nrow(interact.edges)*pcut))
            subpindex = 1:count_ints
            subcounts[i] = median(sorted_ints[subpindex])

            corr_top = [cor(interact.edges[subpindex, p1], interact.edges[subpindex, p2]) for (p1, p2) in combinations(interact.replicate_ids, 2)]
            corr_mean[i] = mean(corr)
            corr_sd[i] = std(corr)
            counts[i] = median(interact.edges.nb_ints[pindex])
            corr_top_mean[i] = mean(corr_top)
        end

        corr = [cor(interact.edges[!, p1], interact.edges[!, p2]) for (p1, p2) in combinations(interact.replicate_ids, 2)]
        corr_mean[length(pcuts)+1] = mean(corr)
        corr_sd[length(pcuts)+1] = std(corr)
        counts[length(pcuts)+1] = median(interact.edges.nb_ints)
        max_count = max(max_count, maximum(counts))
        scatter!(ax_cor, 1:(length(pcuts)+1), corr_mean, label=l, color=color)
        scatter!(ax_count, 1:(length(pcuts)+1), counts, label=l, color=color)
        scatter!(ax_top, 1:length(pcuts), corr_top_mean, label=l, color=color)
        scatter!(ax_ints_count, 1:length(pcuts), subcounts, color=color)
    end
    Legend(gc[1,4], ax_cor, "seed | score")
    Legend(gb[1,4], ax_cor, "seed | score")

    ax2 = Axis(ga[1,3], yticks=(1:length(scores), string.(scores)), xticks = (1:length(seedlens), string.(seedlens)),
        ylabel="minimum alignment score", xlabel="seed length", title="TPR")
    heatmap!(ax2, tpr, colorrange=clims)
    for i in 1:length(seedlens), j in 1:length(scores)
        text!(ax2, string(round(tpr[i,j], digits=2)), position = (i,j), align = (:center, :center), color = tpr[i,j] < 0.5 ? :white : :black)
    end

    ax3 = Axis(ga[1,4], yticks=(1:length(scores), string.(scores)), xticks = (1:length(seedlens), string.(seedlens)),
        ylabel="minimum alignment score", xlabel="seed length", title="FPR")
    heatmap!(ax3, fpr, colorrange=clims)
    for i in 1:length(seedlens), j in 1:length(scores)
        text!(ax3, string(round(fpr[i,j], digits=3)), position = (i,j), align = (:center, :center), color = fpr[i,j] < 0.5 ? :white : :black)
    end

    Colorbar(ga[1,5], limits=clims)

    ax1 = Axis(gb[1,2], title="false discoveries", xlabel="reads per pair", ylabel="false pairs", yscale=log10, xscale=log10)
    ax_ratio = Axis(gb[1,1], title="recovery rate", ylabel="relative recovered reads", xticks=(1:4, ["$(se) | $(ms)" for (se, ms) in first.(params)]), xticklabelrotation = pi/4)
    compare_df = DataFrame(CSV.File(joinpath(assets_path, "csv_resampled/", "stripped_interactions_hfq_hcd.csv")))[1:200, :]
    pair_dict = Dict((row.name1, row.name2)=>row.nb_ints for row in eachrow(compare_df))
    for (i, (((se, ms), label), color)) in enumerate(zip(params, colors))
        resampled_df = DataFrame(CSV.File(joinpath(assets_path, "csv_resampled/", "interactions_resampled_$(se)_$(ms).csv")))
        found_dict = Dict((row.name1, row.name2) => ((row.name1, row.name2) in keys(pair_dict) ? row.nb_ints : 0.0) for row in eachrow(resampled_df))
        false_interactions = Dict{Int, Int}()
        for row in eachrow(resampled_df)
            (row.name1, row.name2) in keys(pair_dict) && continue
            row.nb_ints in keys(false_interactions) ? false_interactions[row.nb_ints] += 1 : false_interactions[row.nb_ints] = 1
        end
        sorted_false_interactions = sort(collect(false_interactions), by=x->x[1])
        scatter!(ax1, first.(sorted_false_interactions), last.(sorted_false_interactions), color=color)

        vals = [v/pair_dict[k] for (k,v) in found_dict if v != 0.0]
        boxplot!(ax_ratio, repeat([i], length(vals)), vals, color=color)
    end

    save("figure_2.svg", fig)

    save( "figure_s1.svg", fig_si)
end
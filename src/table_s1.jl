function checkinteractions(df::DataFrame, verified_pairs::Vector{Tuple{String,String}}; min_reads=2, max_bp_fdr=1.0, max_fisher_fdr=1.0, check_uppercase=true)
	verified_dict = merge(Dict(pair=>[-1, -1, 2.0] for pair in verified_pairs),
							Dict(reverse(pair)=>[-1, -1, 2.0] for pair in verified_pairs))
    check_uppercase && (verified_dict = Dict((uppercase(p[1]), uppercase(p[2]))=>v for (p,v) in verified_dict))
    for row in eachrow(df)
        key = (row[:name1], row[:name2])
        check_uppercase && (key = (uppercase(key[1]), uppercase(key[2])))
        if key in keys(verified_dict)
            verified_dict[key][1] = row[:in_libs]
            verified_dict[key][2] = row[:nb_ints]
            verified_dict[key][3] = row[:bp_pvalue]
        end
    end

	sorted_keys = vcat([[pair, reverse(pair)] for pair in verified_pairs]...)
    check_uppercase && (sorted_keys = [(uppercase(p[1]), uppercase(p[2])) for p in sorted_keys])
	m = reduce(hcat, [verified_dict[key] for key in sorted_keys])'
	verified_stats = DataFrame(
		name1=String[n[1] for n in verified_pairs],
		name2=String[n[2] for n in verified_pairs],
		libs=max.(m[:,1][1:2:end], m[:,1][2:2:end]),
        count=m[:,2][1:2:end] .+ m[:,2][2:2:end],
        fdr=[isnan(first(p)) ? last(p) : (isnan(last(p)) ? first(p) : min(first(p), last(p))) for p in zip(m[:,3][1:2:end] , m[:,3][2:2:end])]
		)
	return verified_stats
end

function checkinteractions(conditions::Vector{DataFrame}, verified_pairs::Vector{Tuple{String,String}}; min_reads=2, max_bp_fdr=1.0, max_fisher_fdr=1.0, check_uppercase=true)
    verified_stats = DataFrame(name1=String[p[1] for p in verified_pairs], name2=String[p[2] for p in verified_pairs])
    for df in conditions
        verified_stats = innerjoin(verified_stats,
            checkinteractions(df, verified_pairs;  min_reads=min_reads, max_bp_fdr=max_bp_fdr, max_fisher_fdr=max_fisher_fdr, check_uppercase=check_uppercase);
                on=[:name1, :name2], makeunique=true)
    end
    return verified_stats
end

function asdataframe(interactions::Interactions; min_reads=3, max_fisher_fdr=1.0, max_bp_fdr=1.0)
    filter_index = (interactions.edges[!, :nb_ints] .>= min_reads) .& (interactions.edges[!, :fisher_fdr] .<= max_fisher_fdr) .&
                        ((interactions.edges[!, :bp_fdr] .<= max_bp_fdr) .| isnan.(interactions.edges.bp_fdr))
    out_df = interactions.edges[filter_index, :]
    out_df[!, :meanlen1] = Int.(round.(out_df[!, :meanlen1]))
    out_df[!, :meanlen2] = Int.(round.(out_df[!, :meanlen2]))
    out_df[!, :nms1] = round.(out_df[!, :nms1], digits=4)
    out_df[!, :nms2] = round.(out_df[!, :nms2], digits=4)
    out_df[:, :name1] = interactions.nodes[out_df[!,:src], :name]
    out_df[:, :name2] = interactions.nodes[out_df[!,:dst], :name]
    out_df[:, :ref1] = interactions.nodes[out_df[!,:src], :ref]
    out_df[:, :ref2] = interactions.nodes[out_df[!,:dst], :ref]
    out_df[:, :type1] = interactions.nodes[out_df[!,:src], :type]
    out_df[:, :type2] = interactions.nodes[out_df[!,:dst], :type]
    out_df[:, :strand1] = interactions.nodes[out_df[!,:src], :strand]
    out_df[:, :strand2] = interactions.nodes[out_df[!,:dst], :strand]
    out_df[:, :left1] = interactions.nodes[out_df[!,:src], :left]
    out_df[:, :left2] = interactions.nodes[out_df[!,:dst], :left]
    out_df[:, :right1] = interactions.nodes[out_df[!,:src], :right]
    out_df[:, :right2] = interactions.nodes[out_df[!,:dst], :right]
    out_df[:, :in_libs] = sum(eachcol(out_df[!, interactions.replicate_ids] .!= 0))
    out_columns = [:name1, :type1, :ref1, :strand1,:left1, :right1, :name2, :type2, :ref2, :strand2, :left2, :right2, :nb_ints, :nb_multi, :in_libs,
    :fisher_pvalue, :fisher_fdr, :odds_ratio, :bp_pvalue, :bp_fdr, :meanlen1, :nms1, :meanlen2, :nms2]
    return sort!(out_df[!, out_columns], :nb_ints; rev=true)
end

function cdsframecoordinate(p::Int, idx::Int, interact::Interactions)
    cds, left, right = interact.nodes[idx, [:cds, :left, :right]]
    tp = interact.nodes.strand[idx] == '-' ? (cds > 0 ? cds : right)-p+1 : p-(cds > 0 ? cds : left)+1
    return tp - Int(tp <= 0)
end

function interaction_set(interact::Interactions, idx1::Int, idx2::Int, window_len::Int)
    set_pairs = Set{Tuple{Int, Int, Int, Int, Float64}}()
    pair = (idx1, idx2)
    if pair in keys(interact.edgestats)
        strand1::Char = interact.nodes.strand[idx1]
        strand2::Char = interact.nodes.strand[idx2]
        for (coord1, coord2) in keys(interact.edgestats[pair][3])
            pv, al1, ar1, al2, ar2, _ = interact.bpstats[(coord1, coord2)]

            al1 = coord1 + (strand1=='-' ? 1 : -1) * (window_len - al1) - Int(strand1=='+')
            ar1 = coord1 + (strand1=='-' ? 1 : -1) * (window_len - ar1)
            al2 = coord2 + (strand2=='-' ? -1 : 1) * (window_len - al2) - Int(strand2=='-')
            ar2 = coord2 + (strand2=='-' ? -1 : 1) * (window_len - ar2)

            cal1 = cdsframecoordinate(al1, idx1, interact)
            car1 = cdsframecoordinate(ar1, idx1, interact)
            cal2 = cdsframecoordinate(al2, idx2, interact)
            car2 = cdsframecoordinate(ar2, idx2, interact)

            push!(set_pairs, (min(cal1, car1), max(cal1, car1), min(cal2, car2), max(cal2, car2), pv))
        end
    end
    return set_pairs
end

function interaction_set(interact::InteractionsNew, idx1::Int, idx2::Int, window_len::Int)
    set_pairs = Set{Tuple{Int, Int, Int, Int, Float64}}()
    pair = (idx1, idx2)
    if pair in keys(interact.edgestats)
        strand1::Char = interact.nodes.strand[idx1]
        strand2::Char = interact.nodes.strand[idx2]
        for (coord1, coord2) in keys(interact.edgestats[pair][3])
            pv, al1, ar1, al2, ar2, _ = interact.bpstats[(idx1, coord1, idx2, coord2)]

            al1 = coord1 + (strand1=='-' ? 1 : -1) * (window_len - al1) - Int(strand1=='+')
            ar1 = coord1 + (strand1=='-' ? 1 : -1) * (window_len - ar1)
            al2 = coord2 + (strand2=='-' ? -1 : 1) * (window_len - al2) - Int(strand2=='-')
            ar2 = coord2 + (strand2=='-' ? -1 : 1) * (window_len - ar2)

            cal1 = cdsframecoordinate(al1, idx1, interact)
            car1 = cdsframecoordinate(ar1, idx1, interact)
            cal2 = cdsframecoordinate(al2, idx2, interact)
            car2 = cdsframecoordinate(ar2, idx2, interact)

            push!(set_pairs, (min(cal1, car1), max(cal1, car1), min(cal2, car2), max(cal2, car2), pv))
        end
    end
    return set_pairs
end

overlap(l1::Int, r1::Int, l2::Int, r2::Int) = min(1.0, max(0.0, min(r2-l1+1, r1-l2+1) / (r1-l1+1)))

function best_overlap(l1::Int, r1::Int, l2::Int, r2::Int, interaction_set::Set{Tuple{Int,Int,Int,Int,Float64}})
    max_overlap1 = 0.0
    max_overlap2 = 0.0
    max_shared_overlap = 0.0
    max_shared_overlap_pvalue = 1.0
    for (al1, ar1, al2, ar2, pv) in interaction_set
        current_overlap1 = overlap(l1, r1, al1, ar1)
        current_overlap2 = overlap(l2, r2, al2, ar2)
        current_shared_overlap = current_overlap1 * current_overlap2
        if current_shared_overlap >= max_shared_overlap
            max_shared_overlap = current_shared_overlap
            max_shared_overlap_pvalue = min(max_shared_overlap_pvalue, pv)
        end
        max_overlap1 = max(current_overlap1, max_overlap1)
        max_overlap2 = max(current_overlap2, max_overlap2)
    end
    length(interaction_set) == 0 && (max_shared_overlap_pvalue = -1.0)
    return max_shared_overlap, max_shared_overlap_pvalue, max_overlap1, max_overlap2
end

function compare_interaction_sites(interact::Interactions, partner1::Vector{String}, partner2::Vector{String},
        from1::Vector{Int}, to1::Vector{Int}, from2::Vector{Int}, to2::Vector{Int}, window_len::Int)
    overlaps = zeros(length(partner1), 4)
    for (i, (n1, n2, l1, r1, l2, r2)) in enumerate(zip(partner1, partner2, from1, to1, from2, to2))
        idx1 = findfirst(interact.nodes.name .== n1)
        idx2 = findfirst(interact.nodes.name .== n2)
        direction1 = interaction_set(interact, idx1, idx2, window_len)
        direction2 = interaction_set(interact, idx2, idx1, window_len)
        olps_1, olps_1_pv, olp1_1, olp2_1 = best_overlap(l1, r1, l2, r2, direction1)
        olps_2, olps_2_pv, olp1_2, olp2_2 = best_overlap(l2, r2, l1, r1, direction2)
        if olps_1 > olps_2
            overlaps[i, 1] = olps_1
            overlaps[i, 2] = olps_1_pv
        else
            overlaps[i, 1] = olps_2
            overlaps[i, 2] = olps_2_pv
        end
        overlaps[i, 3] = max(olp1_1, olp1_2)
        overlaps[i, 4] = max(olp2_1, olp2_2)
    end
    return overlaps
end

function string_class_and_pvalue(overlaps::Matrix{Float64})

    classes = Vector{String}(undef, size(overlaps)[1])
    pvalues = Vector{String}(undef, size(overlaps)[1])
    for (i, r) in enumerate(eachrow(overlaps))
        if r[2] == -1
            classes[i] = "not detected"
            pvalues[i] = "-"
        elseif r[1] > 0.8
            classes[i] = "match"
            pvalues[i] = "$(round(r[2], digits=7))"
        elseif r[1] > 0.0
            classes[i] = "partial"
            pvalues[i] = "$(round(r[2], digits=7))"
        elseif r[3] > 0.8
            classes[i] = "target shifted"
            pvalues[i] = "-"
        elseif r[4] > 0.8
            classes[i] = "sRNA shifted"
            pvalues[i] = "-"
        else
            classes[i] = "no overlap"
            pvalues[i] = "-"
        end
    end
    return classes, pvalues
end

function make_table_s1(assets_folder::String, interact_lcd::Interactions, interact_hcd::Interactions)
    
    df_literature = DataFrame(CSV.File(joinpath(assets_folder, "all_literature.csv"), stringtype=String))
    df_out = DataFrame(sRNA=df_literature.sRNA, target=df_literature.target)
    
    overlaps_hcd = compare_interaction_sites(interact_hcd, df_literature.sRNA, df_literature.target,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_hcd, pvalues_hcd = string_class_and_pvalue(overlaps_hcd)
    df_out.hcd_pvalue = pvalues_hcd
    df_out.hcd_class = classes_hcd
    
    overlaps_lcd = compare_interaction_sites(interact_lcd, df_literature.sRNA, df_literature.target,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_lcd, pvalues_lcd = string_class_and_pvalue(overlaps_lcd)
    df_out.lcd_pvalue = pvalues_lcd
    df_out.lcd_class = classes_lcd
    CSV.write("table_s1.csv", df_out)
end

function make_table_s1_coli(assets_folder::String, interact_il::Interactions, 
        interact_sp::Interactions, interact_lp::Interactions, window_len::Int)

    df_literature = DataFrame(CSV.File(joinpath(assets_folder, "all_literature_coli.csv"), stringtype=String))
    df_out = DataFrame(sRNA=df_literature.sRNA, target=df_literature.target)

    overlaps_il = compare_interaction_sites(interact_il, df_literature.sRNA, df_literature.target,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_il, pvalues_il = string_class_and_pvalue(overlaps_il)
    df_out.ironlimit_pvalue = pvalues_il
    df_out.ironlimit_class = classes_il

    overlaps_sp = compare_interaction_sites(interact_sp, df_literature.sRNA, df_literature.target,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_sp, pvalues_sp = string_class_and_pvalue(overlaps_sp)
    df_out.stationary_pvalue = pvalues_sp
    df_out.stationary_class = classes_sp

    overlaps_lp = compare_interaction_sites(interact_lp, df_literature.sRNA, df_literature.target,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_lp, pvalues_lp = string_class_and_pvalue(overlaps_lp)
    df_out.exponential_pvalue = pvalues_lp
    df_out.exponential_class = classes_lp

    CSV.write("table_s1_coli.csv", df_out)
end
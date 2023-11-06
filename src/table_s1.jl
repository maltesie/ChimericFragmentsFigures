function interaction_counts(interact::Union{Interactions, InteractionsNew}, partners1::Vector{String}, partners2::Vector{String})
    counts = zeros(Int, length(partners1))
    combined_pvalues = ones(Float64, length(partners1)) .* 2.0
    for (i, (n1, n2)) in enumerate(zip(partners1, partners2))
        idx1 = findfirst(interact.nodes.name .== n1)
        idx2 = findfirst(interact.nodes.name .== n2)
        for index_pair in ((idx1, idx2), (idx2, idx1))
            edge_index = findfirst((interact.edges.src .== first(index_pair)) .& (interact.edges.dst .== last(index_pair)))
            if !isnothing(edge_index)
                counts[i] += interact.edges.nb_ints[edge_index]
                !isnan(interact.edges.bp_pvalue[edge_index]) &&
                    (combined_pvalues[i] = min(combined_pvalues[i], interact.edges.bp_pvalue[edge_index]))
            end
        end
    end
	return counts, [pv == 2.0 ? "-" : "$(round(pv, digits=7))" for pv in combined_pvalues]
end

function cdsframecoordinate(p::Int, idx::Int, interact::Union{Interactions, InteractionsNew})
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
            (idx1, coord1, idx2, coord2) in keys(interact.bpstats) || continue #println(interact.nodes.name[idx1], "\n", interact.nodes.name[idx2])
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
    isqrr4 = false
    for (al1, ar1, al2, ar2, pv) in interaction_set
        current_overlap1 = overlap(l1, r1, al1, ar1)
        current_overlap2 = overlap(l2, r2, al2, ar2)
        current_shared_overlap = current_overlap1 * current_overlap2
        #if (al1 == 4) && (ar1 == 27) && (al2 == -144) && (ar2 == -120)
        #    isqrr4 = true
        #    println("current max_overlap: $max_shared_overlap, current_overlap: $current_shared_overlap")
        #end
        if current_shared_overlap >= max_shared_overlap
            max_shared_overlap = current_shared_overlap
            max_shared_overlap_pvalue = (current_shared_overlap == max_shared_overlap) ? min(max_shared_overlap_pvalue, pv) : pv
        end
        max_overlap1 = max(current_overlap1, max_overlap1)
        max_overlap2 = max(current_overlap2, max_overlap2)
    end
    length(interaction_set) == 0 && (max_shared_overlap_pvalue = -1.0)
    #isqrr4 && println(max_shared_overlap)
    return max_shared_overlap, max_shared_overlap_pvalue, max_overlap1, max_overlap2
end

function compare_interaction_sites(interact::Union{Interactions, InteractionsNew}, partner1::Vector{String}, partner2::Vector{String},
        from1::Vector{Int}, to1::Vector{Int}, from2::Vector{Int}, to2::Vector{Int}, window_len::Int)
    overlaps = zeros(length(partner1), 4)
    for (i, (n1, n2, l1, r1, l2, r2)) in enumerate(zip(partner1, partner2, from1, to1, from2, to2))
        idx1 = findfirst(interact.nodes.name .== n1)
        idx2 = findfirst(interact.nodes.name .== n2)
        #(isnothing(idx1) || isnothing(idx2)) && continue
        direction1 = interaction_set(interact, idx1, idx2, window_len)
        direction2 = interaction_set(interact, idx2, idx1, window_len)

        olps_1, olps_1_pv, olp1_1, olp2_1 = best_overlap(l1, r1, l2, r2, direction1)
        olps_2, olps_2_pv, olp1_2, olp2_2 = best_overlap(l2, r2, l1, r1, direction2)
        #if (n1 == "Qrr1") && (n2 == "aphA")
        #    println(i)
        #    println(olps_1)
        #    println(olps_2)
        #    println(direction1)
        #    println(direction2)
        #end
        if olps_1 > olps_2
            overlaps[i, 1] = olps_1
            overlaps[i, 2] = olps_1_pv
        elseif olps_1 == olps_2
            overlaps[i, 1] = olps_1
            overlaps[i, 2] = (olps_1_pv == -1.0) ? olps_2_pv : (olps_2_pv == -1.0 ? olps_1_pv : min(olps_1_pv, olps_2_pv))
        else
            overlaps[i, 1] = olps_2
            overlaps[i, 2] = olps_2_pv
        end
        overlaps[i, 3] = max(olp1_1, olp1_2)
        overlaps[i, 4] = max(olp2_1, olp2_2)
    end
    #println(overlaps[23, :])
    return overlaps
end

function string_class_and_pvalue(overlaps::Matrix{Float64})

    classes = Vector{String}(undef, size(overlaps)[1])
    pvalues = Vector{String}(undef, size(overlaps)[1])
    for (i, r) in enumerate(eachrow(overlaps))
        if r[2] == -1
            classes[i] = "no ligation"
            pvalues[i] = "-"
        elseif r[1] > 0.5
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

function make_table_s1(assets_folder::String, interact_hcd::InteractionsNew, interact_lcd::InteractionsNew, window_len::Int)

    df_literature = DataFrame(CSV.File(joinpath(assets_folder, "all_literature.csv"), stringtype=String))
    df_out = DataFrame(sRNA=df_literature.sRNA, target=df_literature.target)
    df_out.sRNA_region = ["$from - $to" for (from, to) in zip(df_literature.sRNA_from, df_literature.sRNA_to)]
    df_out.target_region = ["$from - $to" for (from, to) in zip(df_literature.target_from, df_literature.target_to)]

    mod_targets = [startswith(target, "vc") ? uppercase(target) : target for target in df_literature.target]

    overlaps_hcd = compare_interaction_sites(interact_hcd, df_literature.sRNA, mod_targets,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)

    classes_hcd, pvalues_hcd = string_class_and_pvalue(overlaps_hcd)
    counts_hcd, combined_pvalues_hcd = interaction_counts(interact_hcd, df_literature.sRNA, mod_targets)
    df_out.hcd_counts = counts_hcd
    #df_out.hcd_cpvalue = combined_pvalues_hcd
    df_out.hcd_pvalue = pvalues_hcd
    df_out.hcd_class = classes_hcd

    #println(overlaps_hcd[23, :])

    overlaps_lcd = compare_interaction_sites(interact_lcd, df_literature.sRNA, mod_targets,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_lcd, pvalues_lcd = string_class_and_pvalue(overlaps_lcd)
    counts_lcd, combined_pvalues_lcd = interaction_counts(interact_lcd, df_literature.sRNA, mod_targets)
    df_out.lcd_counts = counts_lcd
    #df_out.lcd_cpvalue = combined_pvalues_lcd
    df_out.lcd_pvalue = pvalues_lcd
    df_out.lcd_class = classes_lcd

    CSV.write("table_s1.csv", df_out)
end

function make_table_s1_coli(assets_folder::String, interact_il::Interactions,
        interact_sp::Interactions, interact_lp::Interactions, window_len::Int)

    df_literature = DataFrame(CSV.File(joinpath(assets_folder, "all_literature_coli.csv"), stringtype=String))
    df_out = DataFrame(sRNA=df_literature.sRNA, target=df_literature.target)
    df_out.sRNA_region = ["$from - $to" for (from, to) in zip(df_literature.sRNA_from, df_literature.sRNA_to)]
    df_out.target_region = ["$from - $to" for (from, to) in zip(df_literature.target_from, df_literature.target_to)]

    overlaps_il = compare_interaction_sites(interact_il, df_literature.sRNA, df_literature.target,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_il, pvalues_il = string_class_and_pvalue(overlaps_il)
    counts_il, combined_pvalues_il = interaction_counts(interact_il, df_literature.sRNA, df_literature.target)
    df_out.ironlimit_counts = counts_il
    #df_out.ironlimit_cpvalue = combined_pvalues_il
    df_out.ironlimit_pvalue = pvalues_il
    df_out.ironlimit_class = classes_il

    overlaps_sp = compare_interaction_sites(interact_sp, df_literature.sRNA, df_literature.target,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_sp, pvalues_sp = string_class_and_pvalue(overlaps_sp)
    counts_sp, combined_pvalues_sp = interaction_counts(interact_sp, df_literature.sRNA, df_literature.target)
    df_out.stationary_counts = counts_sp
    #df_out.stationary_cpvalue = combined_pvalues_sp
    df_out.stationary_pvalue = pvalues_sp
    df_out.stationary_class = classes_sp

    overlaps_lp = compare_interaction_sites(interact_lp, df_literature.sRNA, df_literature.target,
        df_literature.sRNA_from, df_literature.sRNA_to, df_literature.target_from, df_literature.target_to,
        window_len)
    classes_lp, pvalues_lp = string_class_and_pvalue(overlaps_lp)
    counts_lp, combined_pvalues_lp = interaction_counts(interact_lp, df_literature.sRNA, df_literature.target)
    df_out.exponential_counts = counts_lp
    #df_out.exponential_cpvalue = combined_pvalues_lp
    df_out.exponential_pvalue = pvalues_lp
    df_out.exponential_class = classes_lp

    CSV.write("table_s1_coli.csv", df_out)
end
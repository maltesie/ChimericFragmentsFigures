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
function make_table_s1(assets_folder::String, interact_lcd::Interactions, interact_hcd::Interactions)
    df_literature = DataFrame(CSV.File(joinpath(assets_folder, "all_literature.csv"), stringtype=String))
    df_hcd = asdataframe(interact_hcd)
    df_lcd = asdataframe(interact_lcd)
    df_out = DataFrame(sRNA=df_literature.sRNA, target=df_literature.target)
    df_out.nb_hcd = zeros(Int, nrow(df_out))
    df_out.nb_lcd = zeros(Int, nrow(df_out))
    df_out.fdr_hcd = zeros(Float64, nrow(df_out))
    df_out.fdr_lcd = zeros(Float64, nrow(df_out))
    checked_df = checkinteractions([df_hcd, df_lcd], collect(zip(df_out.sRNA, df_out.target)))
    checked_df.libs[checked_df.libs .== -1] .= 0
    checked_df.libs_1[checked_df.libs_1 .== -1] .= 0
    checked_df.count[checked_df.count .== -2] .= 0
    checked_df.count_1[checked_df.count_1 .== -2] .= 0
    checked_df.fdr[checked_df.fdr .== 2] .= NaN
    checked_df.fdr_1[checked_df.fdr_1 .== 2] .= NaN
    CSV.write("table_s1.csv", checked_df)
end
function make_table_s1_coli(assets_folder::String, interact_il::Interactions, interact_sp::Interactions, interact_lp::Interactions)
    df_literature = DataFrame(CSV.File(joinpath(assets_folder, "all_literature_coli.csv"), stringtype=String))
    df_il = asdataframe(interact_il)
    df_sp = asdataframe(interact_sp)
    df_lp = asdataframe(interact_lp)
    df_out = DataFrame(sRNA=df_literature.sRNA, target=df_literature.target)
    df_out.nb_il = zeros(Int, nrow(df_out))
    df_out.nb_sp = zeros(Int, nrow(df_out))
    df_out.nb_lp = zeros(Int, nrow(df_out))
    df_out.fdr_il = zeros(Float64, nrow(df_out))
    df_out.fdr_sp = zeros(Float64, nrow(df_out))
    df_out.fdr_lp = zeros(Float64, nrow(df_out))
    checked_df = checkinteractions([df_il, df_sp, df_lp], collect(zip(df_out.sRNA, df_out.target)))
    checked_df.libs[checked_df.libs .== -1] .= 0
    checked_df.libs_1[checked_df.libs_1 .== -1] .= 0
    checked_df.libs_2[checked_df.libs_2 .== -1] .= 0
    checked_df.count[checked_df.count .== -2] .= 0
    checked_df.count_1[checked_df.count_1 .== -2] .= 0
    checked_df.count_2[checked_df.count_2 .== -2] .= 0
    checked_df.fdr[checked_df.fdr .== 2] .= NaN
    checked_df.fdr_1[checked_df.fdr_1 .== 2] .= NaN
    checked_df.fdr_2[checked_df.fdr_2 .== 2] .= NaN
    CSV.write("table_s1_coli.csv", checked_df)
end
using Pkg
Pkg.activate( "PlotEnv")
Pkg.instantiate()

using RNASeqTools, BioSequences, BioAlignments, CairoMakie, BioGenerics, DelimitedFiles, GenomicFeatures
using CSV, JLD2, DataFrames, StatsBase, Combinatorics, MultipleTesting, ColorSchemes, HypothesisTests
import GeometryBasics: Polygon

struct Interactions
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int, Dict{Tuple{Int,Int},Int}, Dict{Tuple{Int,Int},Int}}}
    bpstats::Dict{Tuple{Int,Int}, Tuple{Float64, Int64, Int64, Int64, Int64, Float64, Int64}}
    multichimeras::Dict{Vector{Int}, Int}
    replicate_ids::Vector{Symbol}
    counts::Dict{Symbol,Vector{Int}}
end

Interactions(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => Interactions)) do f
    f["interactions"]
end

struct InteractionsNew
    nodes::DataFrame
    edges::DataFrame
    edgestats::Dict{Tuple{Int,Int}, Tuple{Int, Dict{Tuple{Int,Int},Int}, Dict{Tuple{Int,Int},Int}}}
    bpstats::Dict{Tuple{Int, Int, Int, Int}, Tuple{Float64, Int64, Int64, Int64, Int64, Float64, Int64}}
    multichimeras::Dict{Vector{Int}, Int}
    replicate_ids::Vector{Symbol}
    counts::Dict{Symbol,Vector{Int}}
end

InteractionsNew(filepath::String) = jldopen(filepath,"r"; typemap=Dict("ChimericAnalysis.Interactions" => InteractionsNew)) do f
    f["interactions"]
end

data_folder = joinpath(@__DIR__, "data")
source_folder = joinpath(@__DIR__, "src")
interact_hcd = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "hfq_hcd.jld2"))
interact_lcd = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "hfq_lcd.jld2"))
interact_il = Interactions(joinpath(@__DIR__, "data", "jlds", "iron_limit.jld2"))
interact_sp = Interactions(joinpath(@__DIR__, "data", "jlds", "stationary.jld2"))
interact_lp = Interactions(joinpath(@__DIR__, "data", "jlds", "log_phase.jld2"))
interact_clash_exp = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "clash_exp.jld2"))
interact_clash_trans = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "clash_transition.jld2"))
interact_clash_stat = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "clash_stationary.jld2"))
interact_subtilis_ligrseq = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "ligrseq_combined.jld2"))
interact_pseudomonas_stat = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "pseudomonas_stat.jld2"))
interact_pseudomonas_exp = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "pseudomonas_exp.jld2"))
interact_salmonella = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "salmonella.jld2"))
interact_epec_wt = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "epec_wt.jld2"))
interact_epec_cond = InteractionsNew(joinpath(@__DIR__, "data", "jlds", "epec_T3SS_BFP.jld2"))

include(joinpath(source_folder, "figure_1.jl"))
plot_figure_1(joinpath(data_folder, "figure_1"))

include(joinpath(source_folder, "figure_2_figure_s3.jl"))
plot_figure_2(joinpath(data_folder, "figure_2"), [15,25,40], 15, 40, 0, 1, 100000, [12,14,15,16,18], [14,16,17,18,20], 5000000,
    [(12, 14)=>"FPR=0.0143", (15, 14)=>"FPR=0.0210", (12, 17)=>"FPR=0.0013",  (15, 18)=>"FPR=0.0012"])

include(joinpath(source_folder, "figure_3.jl"))
plot_figure_3(joinpath(data_folder, "figure_3"), interact_lcd, 100000, (30,0), (4,5,0,7,8,3), [(12, 14), (15, 14), (12, 17),  (15, 18)])

include(joinpath(source_folder, "figure_4.jl"))
plot_figure_4(joinpath(data_folder, "figure_4"))

include(joinpath(source_folder, "figure_5.jl"))
plot_figure_5(joinpath(data_folder, "figure_5"), interact_hcd)
#all_srna_aggregation_plots(interact_lcd, interact_hcd)

include(joinpath(source_folder, "figure_6.jl"))
plot_figure_6(joinpath(data_folder, "figure_6"), interact_lcd)

include(joinpath(source_folder, "figure_s1.jl"))
plot_figure_s1(joinpath(data_folder, "figure_s1"))

include(joinpath(source_folder, "figure_s2.jl"))
plot_figure_s2(joinpath(data_folder, "figure_s2"))

include(joinpath(source_folder, "figure_s4.jl"))
plot_figure_s4([interact_lcd, interact_hcd], [interact_lp, interact_sp, interact_il],
    [interact_clash_exp, interact_clash_trans, interact_clash_stat], [interact_subtilis_ligrseq], [interact_salmonella],
    [interact_pseudomonas_exp, interact_pseudomonas_stat], [interact_epec_wt, interact_epec_cond],
    [[2400, 1800], [1027, 1844, 1947], [498, 1066, 706], [10462], [1170], [997, 702], [744, 971]], 0.25)

include(joinpath(source_folder, "figure_s5.jl"))
plot_figure_s5(interact_hcd, interact_lcd; bp_interval_len=30, fdr_cut=0.25,
    probmass_min=0.4, probmass_max=0.8, min_partners=3, min_ligpoints=10, connect_within=1)

include(joinpath(source_folder, "figure_s6.jl"))
plot_figure_s6(joinpath(data_folder, "figure_s6"), interact_lcd)

include(joinpath(source_folder, "table_s1.jl"))
make_table_s1(joinpath(data_folder, "table_s1"), interact_hcd, interact_lcd, 30)
make_table_s1_coli(joinpath(data_folder, "table_s1"), interact_il, interact_sp, interact_lp, 30)
make_table_s1_clash(joinpath(data_folder, "table_s1"), interact_clash_exp, interact_clash_trans, interact_clash_stat, 30)

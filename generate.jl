using Pkg
Pkg.activate( "PlotEnv")
Pkg.instantiate()

using RNASeqTools, BioSequences, BioAlignments, CairoMakie, BioGenerics, DelimitedFiles, GenomicFeatures
using CSV, JLD2, DataFrames, StatsBase, Combinatorics, MultipleTesting, ColorSchemes
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

data_folder = joinpath(@__DIR__, "data")
source_folder = joinpath(@__DIR__, "src")
interact_hcd = Interactions(joinpath(@__DIR__, "data", "jlds", "hfq_hcd_12_17.jld2"))
interact_lcd = Interactions(joinpath(@__DIR__, "data", "jlds", "hfq_lcd_12_17.jld2"))


include(joinpath(source_folder, "figure_1.jl"))
plot_figure_1(joinpath(data_folder, "figure_1"))

include(joinpath(source_folder, "figure_2_figure_s1.jl"))
plot_figure_2(joinpath(data_folder, "figure_2"), [15,25,40], 15, 40, 0, 1, 100000, [12,14,15,16,18], [14,16,17,18,20], 2000,
    [(12, 14)=>"FPR=0.0143", (15, 14)=>"FPR=0.0210", (12, 17)=>"FPR=0.0013",  (15, 18)=>"FPR=0.0012"])

include(joinpath(source_folder, "figure_3.jl"))
plot_figure_3(joinpath(data_folder, "figure_3"), interact_lcd, 100000, (30,0), (4,5,0,7,8,3), [(12, 14), (15, 14), (12, 17),  (15, 18)])

include(joinpath(source_folder, "figure_4.jl"))
plot_figure_4(joinpath(data_folder, "figure_4"), interact_lcd)

include(joinpath(source_folder, "figure_5.jl"))
plot_figure_5(joinpath(data_folder, "figure_5"), interact_lcd)

include(joinpath(source_folder, "figure_s2.jl"))
plot_figure_s2(joinpath(data_folder, "figure_s2"))

include(joinpath(source_folder, "figure_s3.jl"))
plot_figure_s3(joinpath(data_folder, "figure_s3"), interact_hcd)

include(joinpath(source_folder, "figure_s4.jl"))
plot_figure_s4(joinpath(data_folder, "figure_s4"), interact_lcd)

include(joinpath(source_folder, "figure_s5.jl"))
plot_figure_s5(joinpath(data_folder, "figure_s5"))

include(joinpath(source_folder, "figure_s6.jl"))
plot_figure_s6(joinpath(data_folder, "figure_s6"))

include(joinpath(source_folder, "table_s1.jl"))
make_table_s1(joinpath(data_folder, "table_s1"), interact_hcd, interact_lcd)
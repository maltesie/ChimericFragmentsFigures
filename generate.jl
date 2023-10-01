using Pkg
Pkg.activate( "PlotEnv")
Pkg.instantiate()

using RNASeqTools, BioSequences, BioAlignments, CairoMakie, BioGenerics, DelimitedFiles, GeometryTypes, CSV, JLD2, DataFrames, StatsBase, Combinatorics

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

data_folder = joinpath(@__DIR__, "data")
source_folder = joinpath(@__DIR__, "src")

include(joinpath(source_folder, "figure_1.jl"))
plot_figure_1(joinpath(data_folder, "figure_1"))

include(joinpath(source_folder, "figure_2_figure_s1.jl"))
plot_figure_2(joinpath(data_folder, "figure_2"), [15,25,40], 15, 40, 0, 1, 100000, [12,14,15,16,18], [14,16,17,18,20], 2000,
    [(12, 14)=>"FPR=0.0143", (15, 14)=>"FPR=0.0210", (12, 17)=>"FPR=0.0013",  (15, 18)=>"FPR=0.0012"])
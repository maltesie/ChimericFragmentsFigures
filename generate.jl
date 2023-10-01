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
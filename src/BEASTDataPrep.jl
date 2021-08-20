module BEASTDataPrep

using CSV, DataFrames, Gadfly, Statistics, PhyloNetworks
using HypothesisTests


export parse_continuous_data,
       plot_transformed_data,
       log_transform_data!,
       standardize_data!,
       conform_tree_and_data,
       conform_traitdata,
       merge_traitdata

include("parse_data.jl")
include("transform_and_scale.jl")
include("trim_taxa.jl")


end

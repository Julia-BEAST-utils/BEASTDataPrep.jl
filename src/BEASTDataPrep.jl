module BEASTDataPrep

using CSV, DataFrames, Gadfly


export parse_continuous_data,
       plot_transformed_data

include("parse_data.jl")
include("transform_and_scale.jl")


end

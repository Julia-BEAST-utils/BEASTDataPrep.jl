module BEASTDataPrep

using CSV, DataFrames, Gadfly, Statistics


export parse_continuous_data,
       plot_transformed_data,
       log_transform_data!,
       standardize_data!

include("parse_data.jl")
include("transform_and_scale.jl")


end

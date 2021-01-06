struct PrepDataFrame
    taxa::Vector{<:AbstractString}
    df::DataFrame

    function PrepDataFrame(df::DataFrame)
        if names(df)[1] != "taxon"
            error("first column of dataframe must be \"$(TAXON_NAME)\"")
        end

        for i = 2:size(df, 2)
            col = df[!, i]
            if !(eltype(col) <: Float64)
                error("all columns (except the first) must have elements of type Float64")
            end
        end

        return new(df[!, 1], df[!, 2:end])
    end
end


import Base.getindex

function getindex(pdf::PrepDataFrame, args...)
    return getindex(pdf.df, args...)
end

import Base.setindex!

function setindex!(pdf::PrepDataFrame, args...)
    return setindex!(pdf.df, args...)
end
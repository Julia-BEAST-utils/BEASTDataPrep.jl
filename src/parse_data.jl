const SPECIES_NAMES = ["taxon", "species"]
const FLOAT_TYPES = [Real, Union{Missing, <:Real}]
const TAXON_NAME = "taxon"

function process_column(col::Vector{<:Real})::Vector{Float64}
    n = length(col)
    return [convert(Float64, x) for x in col]
end

function process_column(col::Union{Missing, <:Real})::Vector{Float64}
    n = length(col)
    v = zeros(n)
    for i = 1:n
        if ismissing(col[i])
            v[i] = NaN
        else
            v[i] = convert(Float64, col[i])
        end
    end
    return v
end


function parse_extension(path::AbstractString)
    last_dot = findlast(isequal('.'), path)
    ext = path[(last_dot + 1):end]
    return ext
end

function parse_continuous_data(datapath::AbstractString)::DataFrame
    ext = parse_extension(datapath)
    if ext in keys(FILE_EXTENSIONS)
        f = FILE_EXTENSIONS[ext]
        return f(datapath)
    end
end

function parse_continuous_data_csv(datapath::AbstractString)
    df = CSV.read(datapath)
    processed_df = parse_continuous_df(df)

    return processed_df
end

function parse_continuous_df(df::DataFrame)
    colnames = names(df)
    taxon_inds = findall(x -> lowercase(x) in SPECIES_NAMES, colnames)
    if length(taxon_inds) > 1
        error("not yet implemented")
    elseif length(taxon_inds) == 0
        error("unable to parse input file. cannot find list of taxa.")
    end

    taxa = df[!, taxon_inds[1]]
    data_inds = Int[]

    processed_df = DataFrame()
    processed_df[TAXON_NAME] = taxa

    nms = names(df)

    for i = 1:size(df, 2)
        col = df[!, i]
        for type in FLOAT_TYPES
            if eltype(col) <: type
                processed_df[nms[i]] = process_column(col)
            end
        end
    end

    return processed_df
end




const FILE_EXTENSIONS = Dict{String, Function}(
                            "csv" => parse_continuous_data_csv
                            )
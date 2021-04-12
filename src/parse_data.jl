using UnPack

const SPECIES_NAMES = ["taxon", "species"]
const FLOAT_TYPES = [Real, Union{Missing, <:Real}]
const TAXON_NAME = "taxon"






function process_column(col::AbstractVector{<:Real})::Vector{Float64}
    n = length(col)
    return [convert(Float64, x) for x in col]
end

function process_column(col::AbstractVector{Union{Missing, T}})::Vector{Float64} where T <: Real
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
    df = DataFrame(CSV.File(datapath))
    processed_df = parse_continuous_df(df)

    return processed_df
end

const TAXA = "TAXA"
const CHARACTERS = "CHARACTERS"
const DATATYPE = "DATATYPE"


function contains_taxon(line::AbstractString)
    line = strip(line)
    return !isempty(line) && line != ";"
end

function parse_taxon_values(s::AbstractString)
    comps = split(strip(s))
    taxon = comps[1]
    values = missing_parse.(comps[2:end])
    return (taxon = taxon, values = values)
end

function missing_parse(s::AbstractString)
    if s in ["?"]
        return NaN
    else
        return parse(Float64, s)
    end
end

# function parse_taxa_nexus(lines::Vector{<:AbstractString}, ind::Int)
#     @assert strip(lines[ind]) == nexus_begin(TAXA)



# function parse_continuous_data_nexus(datapath::AbstractString; title::String = "")
#     lines = readlines(data_path)
#     char_start = -1
#     char_begin = nexus_begin(CHARACTERS)
#     taxa_begin = nexus_begin(TAXA)


#     i = 1
#     while i <= length(lines)
#         line = strip(lines[i])
#         if line == char_begin
#             @unpack title, parse_nexus_character(lines, i)
#         else
#             i+= 1
#         end
# end

const NEXUS_COMMENT = r"\[.*?\]"
const NEXUS_BLOCK = r"\s*BEGIN\s+(.*?);\s*(.*?)\s*END;"s
const NEXUS_TAXA = r"\s*DIMENSIONS\s+NTAX\s*=\s*(\d+?)\s*;.*TAXLABELS\s*(.*?)\s*;"s
const NEXUS_TITLE = r"\s*TITLE\s+(.*?)\s*;"
const NEXUS_DIM = r"DIMENSIONS\s+NCHAR\s*=\s*(\d+)\s*;"
const NEXUS_CHARLABELS = r"\s+CHARSTATELABELS\s+(.*?);"s
const NEXUS_MATRIX = r"\s+MATRIX\s+(.*?)\s+;"s

function nexus_format_regex(s::AbstractString)
    space_or_other = "(?:\\s+|\\s+.+\\s+)"
    r = "FORMAT" * space_or_other * s * "\\s*=\\s*(.*?)(?:\\s+|;)"
    return Regex(r)
end

function parse_characters_nexus(s::AbstractString, n_taxa::Int)
    title = ""
    title_match = match(NEXUS_TITLE, s)
    if !isnothing(title_match)
        title = title_match[1]
    end

    dim = parse(Int, match(NEXUS_DIM, s)[1])

    datatype = match(nexus_format_regex(DATATYPE), s)[1]

    labels_string = match(NEXUS_CHARLABELS, s)[1]
    labels_and_inds = strip.(split(labels_string, ','))
    @assert length(labels_and_inds) == dim
    labels = fill("", dim)
    for i = 1:dim
        split_label = split(labels_and_inds[i])
        @assert parse(Int, split_label[1]) == i
        label = strip(join(split_label[2:end], '_'), ['\'','"'])
        labels[i] = label
    end

    taxa = fill("", n_taxa)
    data = fill(NaN, n_taxa, dim)
    data_lines = split(match(NEXUS_MATRIX, s)[1], '\n')
    @assert length(data_lines) == n_taxa

    for i = 1:n_taxa
        @unpack taxon, values = parse_taxon_values(data_lines[i])
        taxa[i] = taxon
        data[i, :] .= values
    end

    df = [DataFrame(taxon = taxa) DataFrame(data, labels)]

    return (title = title, datatype = datatype, data = df)
end


function split_by_pattern(s::AbstractString, r::Regex)
    matches = RegexMatch[]
    between = String[]

    m = match(r, s)
    start = 1
    while !isnothing(m)
        between_inds = start:(m.offset - 1)
        push!(between, s[between_inds])
        push!(matches, m)

        s = @view s[(m.offset + length(m.match)):end]

        m = match(r, s)
    end
    push!(between, s)
    return (matches = matches, between = between)
end

function remove_comments_nexus(s::AbstractString)
    @unpack between = split_by_pattern(s, NEXUS_COMMENT)
    return join(between, '\n')
end

function parse_continuous_data_nexus(path::String)
    s = read(path, String)
    s = remove_comments_nexus(s)

    @unpack matches = split_by_pattern(s, NEXUS_BLOCK)

    taxon_ind = findall(x -> x[1] == TAXA, matches)
    @assert length(taxon_ind) == 1
    taxa = parse_nexus_taxa(matches[taxon_ind[1]][2])
    n_taxa = length(taxa)

    character_blocks = findall(x -> x[1] == CHARACTERS, matches)

    df = DataFrame(taxon = taxa)

    for block in @view matches[character_blocks]
        content = block[2]
        @unpack datatype, title, data = parse_characters_nexus(content, n_taxa)
        @assert df.taxon == data.taxon
        df = [df data[!, 2:end]]
    end

    return df
end

function parse_nexus_taxa(s::AbstractString)
    m = match(NEXUS_TAXA, s)
    n_taxa = parse(Int, m[1])
    taxa = split(m[2])
    if length(taxa) != n_taxa
        error("could not parse nexus taxa")
    end
    return taxa
end



# function parse_first_block(s::AbstractString)


# function parse_taxa_block(lines::Vector{<:AbstractString})


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
                            "csv" => parse_continuous_data_csv,
                            "nex" => parse_continuous_data_nexus
                            )


function conform_tree_and_data(data::DataFrame, newick::AbstractString;
                               strip_chars::Vector{Char} = Char[],
                               verbose::Bool = true)
    taxa_data = data[!, TAXON_NAME]
    tree = readTopology(newick)
    taxa_tree_original = tipLabels(tree)
    taxa_tree = strip.(taxa_tree_original, Ref(strip_chars))

    # remove taxa that are not present in both tree and data from the tree
    not_in_data = setdiff(taxa_tree, taxa_data)
    for taxon in not_in_data
        if verbose
            @warn "removing $taxon from tree"
        end
        i = findfirst(isequal(taxon), taxa_tree)
        deleteleaf!(tree, taxa_tree_original[i])
    end

    in_both = intersect(taxa_data, taxa_tree)
    in_both_inds = [findfirst(isequal(x), taxa_data) for x in in_both]
    trimmed_df = data[in_both_inds, :]

    return (data = trimmed_df, newick = writeTopology(tree))
end

function merge_traitdata(data::DataFrame...; cut_taxa::Bool = false)
    n = length(data)
    taxon_names = [check_valid(d) for d in data]

    taxon_name = taxon_names[1]
    for i = 2:n
        rename!(data[i], taxon_names[i] => taxon_name)
    end

    args = (on = taxon_name, makeunique = true)

    big_df = cut_taxa ? innerjoin(data...; args...) : outerjoin(data...; args...)
    return big_df
end

function conform_traitdata(data::DataFrame...;kw_args...)
    n = length(data)
    dims = size.(data, 2)

    joint_df = merge_traitdata(data...; kw_args...)
    joint_names = names(joint_df)
    taxon_name = check_valid(joint_df)
    taxa = joint_df[!, taxon_name]

    new_data = Vector{DataFrame}(undef, n)
    ind = 1
    for i = 1:n
        df = DataFrame()
        df[!, taxon_name] = taxa
        for j = 1:(dims[i] - 1) # each of the dims includes a single taxon column
            ind += 1
            df[!, joint_names[ind]] = joint_df[!, ind]
        end

        new_data[i] = df
    end

    return new_data
end

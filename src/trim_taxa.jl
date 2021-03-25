function conform_tree_and_data(data::DataFrame, newick::AbstractString)
    taxa_data = data[!, TAXON_NAME]
    tree = readTopology(newick)
    taxa_tree = tipLabels(tree)

    # remove taxa that are not present in both tree and data from the tree
    not_in_data = setdiff(taxa_tree, taxa_data)
    for taxon in not_in_data
        @warn "removing $taxon from tree"
        deleteleaf!(tree, taxon)
    end

    in_both = intersect(taxa_data, taxa_tree)
    in_both_inds = [findfirst(isequal(x), taxa_data) for x in in_both]
    trimmed_df = data[in_both_inds, :]

    return trimmed_df, writeTopology(tree)
end



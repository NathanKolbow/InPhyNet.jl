# using Distributed


# THINGS TO TEST:
# - all error/warning triggers found below
# - normal runs


function inphynet(
    estgts::AbstractVector{HybridNetwork},
    subsets::AbstractVector{<:AbstractVector{<:AbstractString}};
    verbose::Bool=true,
    snaqargs...
)

    # Gather all unique taxa in the estimated gene trees
    taxa_in_trees = Set{AbstractString}()
    for gt in estgts
        for taxa in tiplabels(gt)
            push!(taxa_in_trees, taxa)
        end
    end
    taxa_in_subsets = Set{AbstractString}(reduce(vcat, subsets))

    # Check inputs for warnings and errors
    if length(taxa_in_trees) != length(taxa_in_subsets)
        taxa_missing_from_subsets = findall(tree_taxa -> !(tree_taxa in taxa_in_subsets), taxa_in_trees)
        taxa_missing_from_trees = findall(subset_taxa -> !(subset_taxa in taxa_in_trees), taxa_in_subsets)
        
        length(taxa_missing_from_subsets) == 0 || @warn "The following taxa were found in `estgts` but not in `subsets`: $(taxa_missing_from_subsets)"
        length(taxa_missing_from_trees) == 0 || throw(ArgumentError("The following taxa were found in `subsets` but do not appear in `estgts`: $(taxa_missing_from_trees)"))
    end
    length(taxa_in_subsets) < sum([length(set) for set in subsets]) && @warn "Some taxa appear in `subsets` multiple times. This may lead to conflicts."
    has_missing_pairs(estgts) && throw(ArgumentError("All pairs of taxa must appear together in at least one gene tree, otherwise pairwise distance is undefined."))
    minimum(length(S) for S in subsets) > 4 || error("All subsets must have length > 4 (required by SNaQ).")


    constraint_networks::AbstractVector{HybridNetwork} = Distributed.pmap(1:length(subsets)) do i
        subset_taxa = subsets[i]

        verbose && println("Beginning pipeline with taxa: $(subset)")

        # 1. Trim gene trees to relevant taxa
        verbose && println("Trimming gene trees to relevant taxa")
        subset_gts = Vector{HybridNetwork}()
        for gt in estgts
            if length(findall(taxa -> taxa in gt)) >= 4
                push!(subset_gts, prune_network(gt, subset_taxa))
            end
        end
        if length(subset_gts) == 0
            throw(ErrorException("0 gene trees contained at least 4 taxa in the following subset: $(subset_taxa)."))
        end

        # 2. subset_gts --> quartet info
        verbose && println("Computing quartet CFs")
        q, t = countquartetsintrees(subset_gts, displayprogress = false)
        df = readTableCF(writeTableCF(q, t))

        # 3. distance matrix
        verbose && println("Computing distance matrix")
        D, namelist =  calculateAGID(subset_gts)

        # 4. snaq starting topology
        verbose && println("Computing NJ network for SNaQ starting topology")
        tre0 = inphynet(D, Vector{HybridNetwork}([]), namelist) # PhyloNetworks.nj! is very, very slow with many taxa

        # 5. SNaQ
        snaq_filename = get(snaqargs, :filename, "")
        if typeof(snaq_filename) <: AbstractString && snaq_filename != ""
            snaq_filename = "$(snaq_filename)$(i)"
        end
        snaq_net = snaq!(tre0, df; snaqargs..., filename = snaq_filename)

        return snaq_net
    end

    verbose && println("Computing distance matrix with full gene trees")
    D, namelist = calculateAGID(estgts)
    return inphynet(D, constraint_networks, namelist), constraint_networks

end


function has_missing_pairs(gts::AbstractVector{HybridNetwork})
    all_taxa = reduce(vcat, [tiplabels(gt) for gt in gts])
    pairs = [(all_taxa[i], all_taxa[j]) for i = 1:(length(all_taxa)-1) for j = (i+1):length(all_taxa)]

    for gt in gts
        gt_taxa = tiplabels(gt)
        
        for j = length(pairs):-1:1
            if pairs[j][1] in gt_taxa && pairs[j][2] in gt_taxa
                deleteat!(pairs, j)
            end
        end

        if length(pairs) == 0 return false end
    end
    return true
end
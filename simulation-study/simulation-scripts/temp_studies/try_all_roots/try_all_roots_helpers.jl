function try_all_roots_netnj!(D::Matrix{Float64}, constraints::Vector{HybridNetwork}, namelist::AbstractVector{<:AbstractString}; 
    supressunsampledwarning=false, copy_retic_names::Bool=false, major_tree_only::Bool=false, force_unrooted::Bool=false)
    
    D = deepcopy(D)
    constraints = deepcopy(constraints)
    namelist = deepcopy(namelist)

    if length(constraints) != 2
        error("ONLY 2 CONSTRAINTS AT A TIME FOR NOW")
        throw(ErrorException())
    end

    PhyloNetworks.check_distance_matrix(D)
    InPhyNet.check_constraints(constraints)
    if !force_unrooted
        InPhyNet.root_constraints!(constraints)
    end
    n = size(D, 1)

    # If namelist is not provided
    if isempty(namelist)
        namelist = string.(1:n)
    end
    if length(namelist) != n
        m = length(namelist)
        throw(ArgumentError("D has dimensions $n x $n but only $m names were provided."))
    elseif length(unique(namelist)) != length(namelist)
        throw(ArgumentError("Names must be unique."))
    end

    # Used for path-finding
    # TODO: `Graph` is a limiting factor in algo speed. If it seems that
    #       we still need to improve algo speed, we can re-use graphs
    #       that are created here instead of making them on-demand
    # full_netgraphs = [Graph(c) for c in constraints]

    # Empty network
    subnets = Vector{InPhyNet.SubNet}([SubNet(i, namelist[i]) for i in 1:n])
    reticmap = InPhyNet.ReticMap(constraints)

    # Edge case: remove (but log) retics that emerge from the root of constraints
    rootretics = Array{Union{Tuple, Edge, Nothing}}(undef, length(constraints))
    for (i, c) in enumerate(constraints)
        hybridbools = [edge.hybrid for edge in c.node[c.root].edge]
        if any(hybridbools)
            sum(hybridbools) == 1 || error("Two reticulations coming out of root, this is not accounted for!")
            hybedge = c.node[c.root].edge[hybridbools][1]
            rootretics[i] = hybedge

            if !supressunsampledwarning
                @warn "Hybridization involving upsampled taxa detected in constraint network $(i). "*
                    "Merging may not behave as expected, see this post for important details: "*
                    "POST NOT MADE YET - PLEASE POST AN ISSUE ON GITHUB"
            end
        else
            # None of the nodes coming out of the root are hybrids, but now
            # we need to make sure that none of the nodes coming out of the
            # root lead to _only_ hybrids, thus effectively making the root
            # lead directly to hybrids.
            children = getchildren(c.node[c.root])
            hybridbools = [[edge.hybrid for edge in child.edge] for child in children]
            needtowarn = supressunsampledwarning

            if sum(hybridbools[1]) == 2
                rootretics[i] = (children[1].edge[1], children[1].edge[2])
            elseif sum(hybridbools[2]) == 2
                rootretics[i] = (children[2].edge[1], children[2].edge[2])
            else
                rootretics[i] = nothing
                needtowarn = false
            end
        end
    end
    rootreticprocessed = [false for _ in 1:length(constraints)]

    # Main algorithm loop
    @debug "----- ENTERING MAIN ALGO LOOP -----"
    while n > 1
        # DEBUG STATEMENT
        # @debug n
        
        possible_siblings = findvalidpairs(constraints, namelist, major_tree_only = major_tree_only, force_unrooted = force_unrooted)
        
        # Find optimal (i, j) idx pair for matrix Q
        i = j = nothing
        try
            i, j = findoptQidx(D, possible_siblings, namelist=namelist) 
        catch e
            if !(typeof(e) <: SolutionDNEError)
                rethrow(e)
            end
        end

        @show n
        i, j = findoptQidx_search(D, constraints, namelist, major_tree_only = major_tree_only, force_unrooted = force_unrooted)
        # i, j = findoptQidx(D, possible_siblings, namelist=namelist)

        if i === nothing
            @debug "Looped through every edge to no avail."
            throw(ErrorException())
        end

        # @debug (namelist[i], namelist[j])
        # @debug "(i, j) = ($(i), $(j))"

        # connect subnets i and j
        subnets[i], edgei, edgej = InPhyNet.mergesubnets!(subnets[i], subnets[j])
        InPhyNet.updateconstraints!(namelist[i], namelist[j], constraints, reticmap, edgei, edgej)

        # if a constraint with a root-retic is down to a single taxa, place that root-retic in the appropriate subnet
        for (cidx, c) in enumerate(constraints)
            if length(c.leaf) == 1 && rootretics[cidx] !== nothing && !rootreticprocessed[cidx]
                fakesubnet = SubNet(-cidx, "__placeholder_for_rootretic_num$(cidx)__")
                subnets[i], edgei, edgej = InPhyNet.mergesubnets!(subnets[i], fakesubnet)
                
                if typeof(rootretics[cidx]) <: Edge
                    InPhyNet.trylogretic!(reticmap, rootretics[cidx], edgei, "from")
                else # type is Tuple
                    InPhyNet.trylogretic!(reticmap, rootretics[cidx][1], edgei, "from")
                    InPhyNet.trylogretic!(reticmap, rootretics[cidx][2], edgej, "from")
                end
                # println("logging retic for edge number $(rootretics[cidx].number)")
                rootreticprocessed[cidx] = true
            end
        end

        # collapse taxa i into j
        for l in 1:n
            if l != i && l != j
                D[l, i] = D[i, l] = (D[l, i] + D[j, l] - D[i, j]) / 2
            end
        end

        # Remove data elements that corresponded to `j`
        idxfilter = [1:(j-1); (j+1):n]
        D = view(D, idxfilter, idxfilter)   # remove j from D
        subnets = view(subnets, idxfilter)
        namelist = view(namelist, idxfilter)

        n -= 1
    end

    mnet = InPhyNet.HybridNetwork(subnets[1].nodes, subnets[1].edges)
    mnet.root = mnet.numNodes
    mnet.node[mnet.root].name = "root"

    mnet = InPhyNet.placeretics!(mnet, reticmap, copy_retic_names=copy_retic_names)
    InPhyNet.removeplaceholdernames!(mnet)

    return mnet
end


function findoptQidx_search(D, constraints, namelist; major_tree_only=false, force_unrooted=false)
    @debug "Looping"
    for (ii, edgei) in enumerate(constraints[1].edge)
        rootonedge!(constraints[1], edgei)
        for (jj, edgej) in enumerate(constraints[2].edge)
            rootonedge!(constraints[2], edgej)

            possible_siblings = findvalidpairs(constraints, namelist, major_tree_only = major_tree_only, force_unrooted = force_unrooted)
            try
                i, j = findoptQidx(D, possible_siblings, namelist=namelist)
                return i, j
            catch e
                if !(typeof(e) <: SolutionDNEError)
                    rethrow(e)
                end
            end
        end
    end
end
# Author: Nathan Kolbow
# Hosted at https://github.com/NathanKolbow/InPhyNet.jl
#
# Source code for main components of the InPhyNet algorithm.



function inphynet_pairwise(D, constraints, namelist; kwargs...)
    D = deepcopy(D)
    constraints = deepcopy(constraints)
    namelist = deepcopy(namelist)

    if length(constraints) == 1
        return inphynet!(D, constraints, namelist)
    end

    for i = 1:(length(constraints) - 1)
        # TODO: do in order of most similar subsets
        @debug "\n\n----------------------------\n\n"
        constraints = constraints[sortperm([c.numtaxa for c in constraints])]

        cs = constraints[[1, 2]]
        leaf_names = Set(union([leaf.name for leaf in cs[1].leaf], [leaf.name for leaf in cs[2].leaf]))
        idxfilter = findall(n -> n in leaf_names, namelist)
        if length(constraints) == 2
            leaf_names = namelist
            idxfilter = 1:length(namelist)
        end

        @debug "--------------------------------------------------"
        @debug "1: $(cs[1].numtaxa), $(cs[1].numhybrids): $(writenewick(cs[1]))"
        @debug "2: $(cs[2].numtaxa), $(cs[2].numhybrids): $(writenewick(cs[2]))"
        temp = inphynet!(deepcopy(view(D, idxfilter, idxfilter)), deepcopy(cs), view(namelist, idxfilter))
        deleteat!(constraints, 2)
        deleteat!(constraints, 1)
        push!(constraints, temp)
    end

    return constraints[1]
end


"""
Runs the InPhyNet algorithm on the given distance matrix and constraint networks
where the entries in `namelist` correspond to indices in `D`.
"""
function inphynet(D::AbstractMatrix{<:Real}, constraints::AbstractVector{HybridNetwork}, namelist::AbstractVector{<:AbstractString}, use_heuristic::Bool = true; kwargs...)
    try
        return inphynet!(deepcopy(D), Vector{HybridNetwork}([readnewick(writenewick(c)) for c in constraints]), deepcopy(namelist); use_heuristic=use_heuristic, kwargs...)
    catch e
        if !(typeof(e) <: ErrorException) || string(e) != "ErrorException(\"No compatible merge found.\")"
            # We don't know what this error is, re-throw it.
            rethrow(e)
        else
            return inphynet_pairwise(D, constraints, namelist)
        end
    end
end
inphynet(D::AbstractMatrix{<:Real}, namelist::AbstractVector{<:AbstractString}; kwargs...) =
    inphynet(D, Vector{HybridNetwork}(), namelist; kwargs...)


function inphynet!(D::AbstractMatrix{<:Real}, constraints::AbstractVector{HybridNetwork}, namelist::AbstractVector{<:AbstractString}; use_heuristic::Bool = true, kwargs...)
    # TODO: re-implement `netnj` with heuristic checking in a manner such that undoing steps is relatively easy
    #       to do so, a lot of the code in the `netnj` body should be compartmentalized - this function, on its
    #       own, should probably around 10-40 lines of code
    # TODO: keep track of pairs of siblings through the algo instead of re-doing work each time - this is something
    #       that should have to be done once and should be very easy to keep track of afterwards

    @debug "Entering inphynet!"
    orig_namelist = deepcopy(namelist)
    orig_constraints = deepcopy(constraints)
    n = size(D, 1)
    @debug "Checking parameters"
    check_inphynet_parameters(D, constraints, namelist; kwargs...)
    
    # Setup
    @debug "Setting up"
    subnets = Vector{SubNet}([SubNet(i, namelist[i]) for i in 1:n])
    reticmap = ReticMap(constraints)
    gammas = [getparentedgeminor(h).gamma for h in collect(keys(reticmap.map))]
    rootretics, rootreticprocessed = setup_root_retics(constraints; kwargs...)
    compatibility_trees = [majortree(c) for c in constraints]
    constraint_sibling_pairs = [findsiblingpairs(c; kwargs...) for c in compatibility_trees]
    
    # Main iterative loop
    while n > 1
        start_time = time_ns()
        @debug n

        possible_siblings = findvalidpairs(compatibility_trees, constraint_sibling_pairs, namelist)
        i, j = findoptQidx(D, possible_siblings, compatibility_trees, namelist = namelist, use_heuristic = use_heuristic)

        subnets[i], edgei, edgej = mergesubnets!(subnets[i], subnets[j])
    
        update_compat_trees!(namelist[i], namelist[j], compatibility_trees, constraint_sibling_pairs)
        updateconstraints!(namelist[i], namelist[j], constraints, reticmap, edgei, edgej; kwargs...)

        fix_root_retics!(constraints, subnets, rootretics, rootreticprocessed, reticmap, i)

        for l = 1:n
            if l != i && l != j
                D[l, i] = D[i, l] = (D[l, i] + D[j, l] - D[i, j]) / 2
            end
        end

        idxfilter = [1:(j-1); (j+1):n]
        D = D[idxfilter, idxfilter]
        subnets = subnets[idxfilter]
        namelist = namelist[idxfilter]
        n -= 1


        time_elapsed = round((time_ns() - start_time) / 1e9, digits=2)
        @debug "Iteration $(n) took $(time_elapsed)s"
    end

    mnet = HybridNetwork(subnets[1].nodes, subnets[1].edges)
    mnet.rooti = mnet.numnodes
    getroot(mnet).name = "root"

    mnet = placeretics!(mnet, reticmap, gammas; kwargs...)
    removeplaceholdernames!(mnet)

    return mnet
end


function fix_root_retics!(constraints::AbstractVector{HybridNetwork}, subnets::AbstractVector{SubNet}, rootretics::AbstractArray, rootreticprocessed::BitVector, reticmap::ReticMap, i)
    # if a constraint with a root-retic is down to a single taxa, place that root-retic in the appropriate subnet
    for (cidx, c) in enumerate(constraints)
        if length(c.leaf) == 1 && rootretics[cidx] !== nothing && !rootreticprocessed[cidx]
            fakesubnet = SubNet(-cidx, "__placeholder_for_rootretic_num$(cidx)__")
            subnets[i], edgei, edgej = mergesubnets!(subnets[i], fakesubnet)
            
            if typeof(rootretics[cidx]) <: Edge
                trylogretic!(reticmap, rootretics[cidx], edgei, "from")
            else # type is Tuple
                trylogretic!(reticmap, rootretics[cidx][1], edgei, "from")
                trylogretic!(reticmap, rootretics[cidx][2], edgej, "from")
            end
            rootreticprocessed[cidx] = true
        end
    end
end


function setup_root_retics(constraints::AbstractVector{HybridNetwork}; supressunsampledwarning::Bool=false)
    rootretics = Array{Union{Tuple, Edge, Nothing}}(undef, length(constraints))
    for (i, c) in enumerate(constraints)
        hybridbools = [edge.hybrid for edge in getroot(c).edge]
        if any(hybridbools)
            sum(hybridbools) == 1 || error("Two reticulations coming out of root, this is not accounted for!")
            hybedge = getroot(c).edge[hybridbools][1]
            rootretics[i] = hybedge

            if !supressunsampledwarning
                @warn "Hybridization involving upsampled taxa detected in constraint network $(i). "*
                    "Merging may not behave as expected, see this post for important details: "*
                    "POST NOT MADE YET - PLEASE POST AN ISSUE ON GITHUB"
            end
        else
            # If there are 0 hybrids in the entire network, we can skip this work
            if length(c.edge) == 0 || sum(edge.hybrid for edge in c.edge) == 0
                rootretics[i] = nothing
                continue
            end

            # None of the nodes coming out of the root are hybrids, but now
            # we need to make sure that none of the nodes coming out of the
            # root lead to _only_ hybrids, thus effectively making the root
            # lead directly to hybrids.
            children = getchildren(getroot(c))
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
    rootreticprocessed = falses(length(constraints))
    return rootretics, rootreticprocessed
end


function check_inphynet_parameters(D::AbstractMatrix{<:Real}, constraints::AbstractVector{HybridNetwork}, namelist::AbstractVector{<:AbstractString}; skip_constraint_check::Bool=false, kwargs...)

    !isempty(namelist) || throw(ArgumentError("Please provide a non-empty `namelist`."))
    PhyloNetworks.check_distance_matrix(Matrix{Float64}(D))
    if skip_constraint_check
        @warn "Skipping constraint validation - unexpected behavior may occur if one or more constraints violate certain conditions."
    else
        check_constraints!(constraints; kwargs...)
    end

    if length(namelist) != size(D, 1)
        throw(ArgumentError("`namelist` has length $(length(namelist)) but `D` has size ($(size(D, 1)), $(size(D,2)))."))
    end

    if length(unique(namelist)) != length(namelist)
        throw(ArgumentError("Namelist contains non-unique names."))
    end
end


@inline function removeplaceholdernames!(mnet::HybridNetwork)
    for node in mnet.node
        if startswith(node.name, "__placeholder_for_rootretic_num")
            deleteleaf!(mnet, node)
        elseif !node.leaf && startswith(node.name, "___internal")
            node.name = ""
        end
    end
end


"""

Checks validity of input constraints. So far, the only check is to make sure
that all nodes have exactly 3 edges except for the root.
"""
function check_constraints!(constraints::Vector{HybridNetwork}; kwargs...)
    for (i, constraint) in enumerate(constraints)
        check_constraint!(i, constraint, true; kwargs...)
    end
end


"""

If any constraints are unrooted, picks a root node from any node directly
under the root (somewhat randomly).
"""
function root_constraints!(constraints::Vector{HybridNetwork})
    for c in constraints
        root_edges = getroot(c).edge
        if length(root_edges) > 2
            if any(getchild(E).leaf for E in root_edges)
                rootonedge!(c, root_edges[findfirst(E -> getchild(E).leaf, root_edges)])
            else
                for root_edge in root_edges
                    try
                        rootonedge!(c, root_edge)
                        break
                    catch
                    end
                end
            end
        end
    end
end


"""

Checks validity of a single input constraint networks. Checks include:
1. All nodes have exactly 3 edges except the root (unless the network is a single taxa)
2. Reticulations do not lead directly into other reticulations
"""
function check_constraint!(idx::Int64, net::HybridNetwork, autofix::Bool=true; depth::Int=0, kwargs...)
    if net.numtaxa == 1 return end

    # Check #1
    for (node_idx, node) in enumerate(net.node)
        if node.leaf
            if length(node.edge) != 1
                throw(ConstraintError(idx, "Leaf nodes must have exactly 1 attached edge."))
            end
        elseif node == getroot(net)
            if length(node.edge) != 2
                if autofix && depth < 5
                    root_constraints!([net])
                    return check_constraint!(idx, net, true, depth=depth+1)
                else
                    throw(ConstraintError(idx, "Root node must have 2 attached edges (even if net is treated as unrooted - no polytomies allowed)."))
                end
            end
        elseif length(node.edge) != 3
            if autofix && length(node.edge) == 2 && depth < 5
                PhyloNetworks.fuseedgesat!(node_idx, net)
            else
                throw(ConstraintError(idx, "Internal nodes must have exactly 3 attached edges."))
            end
        end
    end

    # Check #2: try to make it so that the root does NOT have a reticulation coming out of it
    i = 0
    edges = net.edge
    retic_at_root = any(e.hybrid for e in getroot(net).edge)
    while retic_at_root
        i += 1

        if i > length(edges)
            throw(ConstraintError(idx, "Could not find a valid root placement where a reticulation is not emerging from the root."))
        end

        try
            rootonedge!(net, edges[i])
            retic_at_root = any(e.hybrid for e in getroot(net).edge)
            if !retic_at_root break end
        catch
        end
    end

end
check_constraint!(net::HybridNetwork; kwargs...) = check_constraint!(0, net; kwargs...)


"""

Updates `net` to include the reticulations that we've kept track of along the way
in our algo but haven't placed yet.
"""
function placeretics!(net::HybridNetwork, reticmap::ReticMap, gammas; copy_retic_names::Bool=false, kwargs...)
    namepairs = []
    retic_names = []
    counter = 0
    
    check_reticmap(reticmap)
    collected_keys = collect(keys(reticmap.map))
    collected_keys = collected_keys[sortperm([k.name for k in collected_keys])]

    for hyb in collected_keys
        from = reticmap.map[hyb][1]
        to = reticmap.map[hyb][2]
        
        if from.node[1].name == ""
            from.node[1].name = "___internal"*string(counter+=1)
        end
        if from.node[2].name == ""
            from.node[2].name = "___internal"*string(counter+=1)
        end
        
        if to.node[1].name == ""
            to.node[1].name = "___internal"*string(counter+=1)
        end
        if to.node[2].name == ""
            to.node[2].name = "___internal"*string(counter+=1)
        end

        # Alphanumeric name ordering for easier uniqueness checking
        name1a, name1b = from.node[1].name < from.node[2].name ? [from.node[1].name, from.node[2].name] : [from.node[2].name, from.node[1].name]
        name2a, name2b = to.node[1].name < to.node[2].name ? [to.node[1].name, to.node[2].name] : [to.node[2].name, to.node[1].name]
        newpair = ([name1a, name1b], [name2a, name2b])

        if !(newpair in namepairs)
            push!(namepairs, newpair)
            push!(retic_names, hyb.name)
        end
    end
    mnet = readnewick(writenewick(net))

    edgepairs = []
    for ((fromname1, fromname2), (toname1, toname2)) in namepairs
        fromedge = nothing
        toedge = nothing

        for edge in mnet.edge
            nodenames = [n.name for n in edge.node]
            if fromname1 in nodenames && fromname2 in nodenames
                fromedge = edge
            elseif toname1 in nodenames && toname2 in nodenames
                toedge = edge
            end
        end
        if fromedge === nothing
            error("from edge nothing")
        elseif toedge === nothing
            error("to edge nothing")
        end
        push!(edgepairs, [fromedge, toedge])
    end

    i = 1
    while i <= length(edgepairs)
        fromedge, toedge = edgepairs[i]

        hybnode, hybedge = addhybridedge!(mnet, fromedge, toedge, true)
        getparentedgeminor(hybnode).gamma = gammas[i]
        getparentedge(hybnode).gamma = 1 - gammas[i]

        hybnode.name = copy_retic_names ? retic_names[i] : hybnode.name
        mnet.rooti = findfirst([n.name == "root" for n in mnet.node])

        for j=(i+1):length(edgepairs)
            if edgepairs[2] == toedge
                edgepairs[2] = hybedge
            end
        end
        i += 1
    end

    return mnet
end


"""

Updates constraint networks after (i, j) with names (nodenamei, nodenamej) have been
merged. Also updates `reticmap` to keep track of any reticulations that get removed
in this process.

By convention we keep `nodenamei` and replace node names with `nodenamej`
"""
function updateconstraints!(nodenamei::AbstractString, nodenamej::AbstractString, 
    constraints::Vector{HybridNetwork}, reticmap::ReticMap,
    subnetedgei::Edge, subnetedgej::Edge; kwargs...)

    for (netidx, net) in enumerate(constraints)
        idxi = -1
        idxj = -1
        for (nodeidx, node) in enumerate(net.node)
            if node.name == nodenamei idxi = nodeidx
            elseif node.name == nodenamej idxj = nodeidx
            elseif idxi != -1 && idxj != -1 break
            end
        end

        hasi = idxi != -1
        hasj = idxj != -1
        if hasj && !hasi
            net.node[idxj].name = nodenamei
            
            # Don't need to update sibling pairs here b/c changing the name does that for us
        elseif hasi && hasj
            nodei = net.node[idxi]
            nodej = net.node[idxj]

            @debug "Merging in net #$(netidx)"
            mergeconstraintnodes!(net, nodei, nodej, reticmap, subnetedgei, subnetedgej)
        end
    end
end


function update_compat_trees!(nodenamei::AbstractString, nodenamej::AbstractString, compat_trees::AbstractVector{HybridNetwork}, constraint_sibling_pairs::AbstractVector)
    
    for (netidx, net) in enumerate(compat_trees)
        idxi = -1
        idxj = -1
        for (nodeidx, node) in enumerate(net.node)
            if node.name == nodenamei idxi = nodeidx
            elseif node.name == nodenamej idxj = nodeidx
            elseif idxi != -1 && idxj != -1 break
            end
        end

        hasi = idxi != -1
        hasj = idxj != -1

        if hasj && !hasi
            net.node[idxj].name = nodenamei
        elseif hasi && hasj
            if net.numtaxa > 2
                deleteleaf!(net, net.node[idxj])
            else
                compat_trees[netidx] = prune_network(net, [nodenamei])
            end

            constraint_sibling_pairs[netidx] = findsiblingpairs(net)
        end
    end
end


function mergeconstraintnodes!(net::HybridNetwork, nodei::Node, nodej::Node, reticmap::ReticMap, subnetedgei::Edge, subnetedgej::Edge)
    parentsi = getnodes(nodei)
    parentsj = getnodes(nodej)

    length(parentsi) == 1 || length(net.node) == 1 || error("Found $(length(parentsi)) nodes above leaf $(nodei.name)?")    # sanity check, remove when finalized
    length(parentsj) == 1 || length(net.node) == 1 || error("Found $(length(parentsj)) nodes above leaf $(nodej.name)?")    # sanity check; remove when finalized
    
    parentsi = parentsi[1]
    parentsj = parentsj[1]

    if (parentsi == parentsj && length(net.leaf) == 2) || (parentsi == nodej && parentsj == nodei)
        @debug "a: ($(nodei.name), $(nodej.name))"

        # TODO: clean this up, when they're nothing we're assigning them randomly right now
        for edge in net.edge
            if edge.hybrid && !edge.ismajor
                try
                    if reticmap.map[getchild(edge)][1] === nothing
                        logretic!(reticmap, edge, subnetedgei, "from")
                    end
                    if reticmap.map[getchild(edge)][2] === nothing
                        logretic!(reticmap, edge, subnetedgej, "to")
                    end
                catch end
            end
        end
        ################

        for edge in net.edge deleteEdge!(net, edge) end
        deleteNode!(net, nodej)
        net.node = [nodei]
        net.edge = []
        net.numtaxa = 1
        net.numnodes = 1
        net.numedges = 0
        net.rooti = 1
        nodei.edge = []
    elseif length(net.leaf) == 2
        @debug "e: ($(nodei.name), $(nodej.name))"

        # Some retics may still exist in the net, but only 2 leaves are left, so clean up
        nodes_above_nodei = []
        nodes_above_nodej = []
        search_node = nodei
        while search_node != getroot(net)
            push!(nodes_above_nodei, search_node)
            try
                search_node = getparent(getparentedge(search_node))
            catch
                search_node = getparent(getparentedgeminor(search_node))
            end
        end
        search_node = nodej
        while search_node != getroot(net)
            push!(nodes_above_nodej, search_node)
            try
                search_node = getparent(getparentedge(search_node))
            catch
                search_node = getparent(getparentedgeminor(search_node))
            end
        end
        push!(nodes_above_nodei, getroot(net))
        push!(nodes_above_nodej, getroot(net))
        mrca = intersect(nodes_above_nodei, nodes_above_nodej)[1]

        # Find all retics that come from nodei
        search_node = nodei
        hybedges_i = []
        hybedges_j = []
        # while !(mrca in getparents(search_node))
        while search_node != mrca
            if length(getparents(search_node)) != 1
                search_node = getparent(getparentedge(search_node))
            else
                try
                    search_node = getparent(search_node)
                catch
                    search_node = getparent(getparentedgeminor(search_node))
                end
            end

            for edge in search_node.edge
                if edge.hybrid # && !edge.ismajor
                    if getchild(edge) == search_node
                        push!(hybedges_i, (edge, "to"))
                    else
                        push!(hybedges_i, (edge, "from"))
                    end
                end
            end
        end
        search_queue = [search_node] # now search going downwards to find all children
        while length(search_queue) > 0
            search_node = pop!(search_queue)
            for edge in search_node.edge
                if edge.hybrid # && !edge.ismajor
                    if !((edge, "from") in hybedges_i)
                        push!(hybedges_i, (edge, "from"))
                    else
                        push!(hybedges_j, (edge, "to"))
                    end
                end
            end

            for child in getchildren(search_node) push!(search_queue, child) end
        end

        # Find all retics that come from nodej
        search_node = nodej
        # while !(mrca in getparents(search_node))
        while search_node != mrca
            if length(getparents(search_node)) != 1
                search_node = getparent(getparentedge(search_node))
            else
                try
                    search_node = getparent(search_node)
                catch
                    search_node = getparent(getparentedgeminor(search_node))
                end
            end

            for edge in search_node.edge
                if edge.hybrid && !edge.ismajor
                    if getchild(edge) == search_node
                        push!(hybedges_j, (edge, "to"))
                    else
                        push!(hybedges_j, (edge, "from"))
                    end
                end
            end
        end
        search_queue = [search_node] # now search going downwards to find all children
        while length(search_queue) > 0
            search_node = pop!(search_queue)
            for edge in search_node.edge
                if edge.hybrid && !edge.ismajor
                    if !((edge, "from") in hybedges_j)
                        push!(hybedges_j, (edge, "from"))
                    else
                        push!(hybedges_i, (edge, "to"))
                    end
                end
            end

            for child in getchildren(search_node) push!(search_queue, child) end
        end
        
        hybedges_i = unique(hybedges_i)
        hybedges_j = unique(hybedges_j)

        # Process `from` edges first
        for (edge, direction) in hybedges_i
            if direction != "from" continue end
            trylogretic_single!(reticmap, edge, subnetedgei, "from")
            # logretic!(reticmap, edge, subnetedgei, "from")
        end
        for (edge, direction) in hybedges_j
            if direction != "from" continue end
            trylogretic_single!(reticmap, edge, subnetedgej, "from")
        end
        
        # Now process `to` edges
        for (edge, direction) in hybedges_i
            if direction != "to" continue end
            trylogretic_single!(reticmap, edge, subnetedgei, "to")
        end
        for (edge, direction) in hybedges_j
            if direction != "to" continue end
            trylogretic_single!(reticmap, edge, subnetedgej, "to")
        end

        #
        for edge in net.edge deleteEdge!(net, edge) end
        deleteNode!(net, nodej)
        net.node = [nodei]
        net.edge = []
        net.numtaxa = 1
        net.numnodes = 1
        net.numedges = 0
        net.rooti = 1
        nodei.edge = []
    elseif parentsi == parentsj && (parentsi == getroot(net) || length(getchildren(parentsi)) > 2)
        # this case happens when `net` is unrooted and nodei & nodej are "outgroup" taxa
        # e.g., this case would happen when joining a & b in (a, b, c, d)
        deleteNode!(net, nodej)
        for edge in net.edge
            deleteEdge!(net, edge)
        end
        # deleteleaf!(net, nodej)
    elseif parentsi == parentsj
        @debug "b: ($(nodei.name), $(nodej.name))"
        # 04/25/2024:   (this note is written well after this code was first written)
        #               this code could likely be simplified almost entirely
        #               to just `deleteleaf!(nodei)`

        # no reticulations: just merge the nodes
        parent = parentsi
        
        # 1. remove the edge connecting nodej and its parent
        parent.edge = filter(e -> e != nodej.edge[1], parent.edge)
        deleteEdge!(net, nodej.edge[1], part=false)
        nodej.edge = []

        # 2. remove node j entirely
        deleteNode!(net, nodej)

        internal = [node for node in getnodes(parent) if node != nodei && node != nodej][1]
        edge1 = nodei.edge[1]
        edge2 = ifelse(parent.edge[1] == edge1, parent.edge[2], parent.edge[1])

        # 3. the links are now: (internal) --edge2-- (parent) --edge1-- (nodei)
        #              we want: (internal) --edge2-- (nodei)
        #    a) remove parent
        #    b) remove edge1
        #    c) remove edge2
        #    c) connect internal and nodei w/ a new edge
        parent.edge = []
        deleteNode!(net, parent)
        
        deleteEdge!(net, edge1)
        deleteEdge!(net, edge2)
        internal.edge = filter(e -> e != edge1 && e != edge2, internal.edge)
        nodei.edge = filter(e -> e != edge1 && e != edge2, nodei.edge)

        newedge = connectnodes!(nodei, internal)  # handy fxn from SubNet.jl
        push!(net.edge, newedge)
    elseif major_mrca(nodei, nodej, getroot(net)) == getroot(net)
        # Merging unrooted constraint across the "root"
        if length(net.node) == 3
            @debug "cross-root B: ($(nodei.name), $(nodej.name))"
            # tricky case b/c ((a, b), c) all pairs of taxa are valid siblings
            # so keep retics straight is tough
            throw(ErrorException("Case not implemented yet."))
        else
            @debug "cross-root A: ($(nodei.name), $(nodej.name))"
            graph, W, nodesinpath, edgesinpath = find_valid_node_path(net, nodei, nodej)

            if any(i -> !isassigned(nodesinpath, i), 1:length(nodesinpath))
                throw(ErrorException("Runtime error: no path exists connecting $(nodei.name) and $(nodej.name)"))
            end

            # Log retics
            relevanttoi = true
            logged_edgeinpath = false
            for (node_idx, node) in enumerate(nodesinpath)
                logged_edgeinpath = log_edge_path_retics_from_node(
                    node, edgesinpath, relevanttoi, false, nothing,
                    reticmap, subnetedgei, subnetedgej,
                    net, logged_edgeinpath, nodesinpath, nodei, nodej
                )

                if node == getroot(net)
                    relevanttoi = false
                end
            end

            # Delete the node that is alone on the other side of the root
            curr = getparent(nodei)
            if has_direct_root_connection(net, nodej)
                curr = getparent(nodej)
                getparent(nodej).edge = filter(e -> e != getparentedge(nodej), getparent(nodej).edge)

                deleteNode!(net, nodej)
                e = getparentedge(nodej)
                deleteEdge!(net, e)
                for node in e.node node.edge = filter(node_edge -> node_edge != e, node.edge) end
            elseif has_direct_root_connection(net, nodei)
                nodei_name = nodei.name
                curr = getparent(nodei)
                getparent(nodei).edge = filter(e -> e != getparentedge(nodei), getparent(nodei).edge)

                deleteNode!(net, nodei)
                e = getparentedge(nodei)
                deleteEdge!(net, e)
                for node in e.node node.edge = filter(node_edge -> node_edge != e, node.edge) end
                nodej.name = nodei_name
            else
                throw(ErrorException("Expected at least one of the node's parents to be the root. Unknown case."))
            end


            # If there are a series of nodes above `nodej` that are not the node (nodes leading to retics): remove those nodes
            while curr != getroot(net) && length(getparents(curr)) > 0
                n_child = length(getchildren(curr))
                if !(n_child == 0 || (any(edge.hybrid for edge in curr.edge) && n_child == 1))
                    break
                end
                @debug "Removing leftover nodes"

                next_curr = getparents(curr)
                if length(next_curr) == 1
                    next_curr = next_curr[1]
                elseif length(next_curr) == 2
                    next_curr = getparent(getparentedge(curr))
                end
                
                getparent(curr).edge = filter(e -> e != getparentedge(curr), getparent(curr).edge)
                deleteNode!(net, curr)
                deleteEdge!(net, getparentedge(curr))
                curr = next_curr
            end

            # Disconnect any hybrid edges that were logged
            for node in nodesinpath
                edges = node.edge
                isminorhyb = [e.hybrid && !e.ismajor for e in edges]

                for hyb_edge in edges[isminorhyb]
                    @debug "Disconnecting edge"
                    for node in hyb_edge.node
                        if !(node in nodesinpath) continue end
                        node.edge = filter(e -> e != hyb_edge, node.edge)
                    end

                    # Changing an edge's node list to only have 1 node causes errors, so replace w/ a placeholder
                    if hyb_edge.node[1] in nodesinpath
                        hyb_edge.node[1] = PhyloNetworks.Node(abs(rand(Int64)), false, false)
                    end
                    if hyb_edge.node[2] in nodesinpath
                        hyb_edge.node[2] = PhyloNetworks.Node(abs(rand(Int64)), false, false)
                    end
                end
            end
        end

        # If the root node is redundant, move the root node down
        while length(getroot(net).edge) == 1 && !getchild(getroot(net)).leaf
            @debug "Moving root node down"
            old_root = getroot(net)
            net.rooti = findfirst(node -> node == getchild(getroot(net)), net.node)
            
            deleteNode!(net, old_root)
            for node in getnodes(old_root)
                node.edge = filter(e -> old_root != e.node[1] && old_root != e.node[2], node.edge)
            end
        end
    else
        @debug "c: ($(nodei.name), $(nodej.name))"
        graph, W, nodesinpath, edgesinpath = find_valid_node_path(net, nodei, nodej)

        # find the node that should be the new tip after merging
        newtip = nodesinpath[1]
        while length(getparents(newtip)) > 0
            parents = getparents(newtip)
            if length(parents) == 2
                if parents[1] in nodesinpath
                    newtip = parents[1]
                    continue
                elseif parents[2] in nodesinpath
                    newtip = parents[2]
                    continue
                else
                    break
                end
            elseif !(parents[1] in nodesinpath)    # length(parents) == 1
                break
            end

            newtip = getparent(newtip)
        end
        
        ishybedge = [e.hybrid && !e.ismajor for e in edgesinpath]
        hybedge = nothing
        hybedgelogged = false
        if newtip === nothing && sum(ishybedge) == 1
            # we follow a hybrid at some point; the new tip is the shared parent of the node at the beginning
            # and the node at the end of the hybrid edge
            hybedge = edgesinpath[ishybedge][1]

            if getparent(hybedge.node[1]) == getparent(hybedge.node[2])
                newtip = getparent(hybedge.node[1])
                for edge in newtip.edge
                    if hybedge.node[1] in edge.node
                        push!(edgesinpath, edge)
                    elseif hybedge.node[2] in edge.node
                        push!(edgesinpath, edge)
                    end
                end
            else
                if length(getparents(hybedge.node[1])) == 2
                    newtip = hybedge.node[2]
                    for e in hybedge.node[1].edge
                        if !(e in edgesinpath)
                            push!(edgesinpath, e)
                        end
                    end
                else
                    newtip = hybedge.node[1]
                    for e in hybedge.node[2].edge
                        if !(e in edgesinpath)
                            push!(edgesinpath, e)
                        end
                    end
                end
            end

            # TODO: direction shouldn't lead to any rooting issues b/c
            # we're at the leaves, but direction not being retained here is bad
            hybedgelogged = true
            logretic!(reticmap, hybedge, subnetedgei, "from")
            logretic!(reticmap, hybedge, subnetedgej, "to")
        elseif newtip === nothing && any(ishybedge)
            error("Unknown case. Multiple hybrid edges in the path merging $(nodei.name) and $(nodej.name).")
        elseif newtip === nothing
            error("Unknown case. No newtip found but we also don't take a hybrid path.")
        end

        # find the retics that we need to keep track of
        # our a_star path goes from `i` to `j`, so any retics we find are relevant to `i` up
        # until we cross the new tip, at which point they become relevant to `j`
        relevanttoi = true
        edgeinpath_logged = false
        for (node_idx, node) in enumerate(nodesinpath)
            # Log all retics at this point
            edgeinpath_logged = log_edge_path_retics_from_node(
                node, edgesinpath, relevanttoi, hybedgelogged, hybedge,
                reticmap, subnetedgei, subnetedgej,
                net, edgeinpath_logged, nodesinpath, nodei, nodej
            )

            if node == newtip
                relevanttoi = false
            end
        end

        ############################
        ######## EDGE CASES ########
        ############################

        # if any nodes are (1) a hybrid, (2) just above a leaf, (3) have a hybrid just above its minor edge parent,
        # they are an edge case
        relevanttoi = true
        for node in nodesinpath
            # Only interested in hybrid nodes
            if !node.hybrid
                if node == newtip relevantoi = false end
                continue
            end

            # Skip if it doesn't have a minor parent
            try getparentedgeminor(node) catch e
                if node == newtip relevantoi = false end
                continue
            end

            buggy_ancestor = getparent(getparentedgeminor(node))
            if any(child.leaf for child in getchildren(node)) && any(parent.hybrid for parent in getparents(buggy_ancestor))
                parent_hybrid = getparents(buggy_ancestor)[findfirst(parent.hybrid for parent in getparents(buggy_ancestor))]
                hyb_edge = nothing
                try
                    # If this fails, then the minor edge has already been removed (i.e. logged) and we can just move on
                    hyb_edge = getparentedgeminor(parent_hybrid)
                catch
                    if node == newtip relevantoi = false end
                    continue
                end
                trylogretic_single!(reticmap, hyb_edge, ifelse(relevanttoi, subnetedgei, subnetedgej), "to")
            end
            if node == newtip relevantoi = false end
        end

        ############################
        ############################
        ############################

        # purge all the nodes we've passed through to this point from the network
        for node in nodesinpath
            if node == newtip continue end
            if node.hybrid && !fully_logged(reticmap, node) continue end
            deleteNode!(net, node)
        end
        for edge in edgesinpath
            deleteEdge!(net, edge)
            for node in edge.node
                node.edge = filter(e -> e != edge, node.edge)
            end
        end

        newtip.leaf = true
        newtip.hybrid = false
        push!(net.leaf, newtip)
        net.numtaxa += 1
        newtip.name = nodei.name

        fuseredundantedges!(net)
    end
end


"""

Gets the MRCA of `nodei` and `nodej` along their major tree.
"""
function major_mrca(nodei::Node, nodej::Node, root::Node)
    obs_nodes = Set()
    while true
        if nodei in obs_nodes && nodei != root
            return nodei
        elseif nodej in obs_nodes && nodej != root
            return nodej
        elseif nodei == nodej
            return nodei
        end

        push!(obs_nodes, nodei)
        push!(obs_nodes, nodej)

        parentsi = getparents(nodei)
        if length(parentsi) == 1
            nodei = parentsi[1]
        elseif length(parentsi) == 2
            nodei = getparent(getparentedge(nodei))
        end

        parentsj = getparents(nodej)
        if length(parentsj) == 1
            nodej = parentsj[1]
        elseif length(parentsj) == 2
            nodej = getparent(getparentedge(nodej))
        end
    end
end


"""

Helper function for when an edge (`notinpath_edge`) appears on but not within a merge path.
Finds and returns:
    1. Whether the edge should be logged with direction "to" or "from"
    2. Whether the edge should be associated with `subnetedgei` or `subnetedgej`
"""
function find_not_in_path_edge_vector(notinpath_edge::Edge, subnetedgei::Edge, subnetedgej::Edge, nodesinpath::Array{Node}, nodei::Node, nodej::Node)
    
    notinpath_direction = (getchild(notinpath_edge) in nodesinpath) ? "to" : "from"

    # Find MRCA of nodei and nodej. The relevant subnet edge depends on which side of the MRCA the notinpath_edge is on.
    relevant_edge = subnetedgej
    mrca = nodei
    while true
        parents = getparents(mrca)
        if !any(parent in nodesinpath for parent in parents)
            break
        end

        # If we find `notinpath_edge`, we found it from `i` so set relevant subnet edge to `i`.
        for edge in mrca.edge
            if edge == notinpath_edge
                relevant_edge = subnetedgei
            end
        end

        # Go to next parent
        for parent in parents
            if parent in nodesinpath mrca = parent end
        end
    end

    return notinpath_direction, relevant_edge
end


"""

Helper function to remove redundant edges (e.g. A --> (internal node) --> (internal node))
that arise from one case of the function `mergeconstraintnodes!`.

Redundant nodes in this case will always have only two edges while *not* being the root.
"""
function fuseredundantedges!(net::HybridNetwork)
    for (i, node) in enumerate(net.node)
        if node.leaf || node.hybrid || node == getroot(net) continue end
        if length(node.edge) == 2
            fuseedgesat!(i, net)
        end
    end
end


"""

Helper function used at one point in `mergeconstraintnodes!`. Gets all of the descendants of
`node` that do not have any children and are not in `path`
"""
function gather_hyb_descendants_outside_of_path(node::Node, path::Array{Node}; stopathyb::Bool=false)
    val = gather_hyb_descendants_outside_of_path_recur(node, path, stopathyb=stopathyb)
    if val === nothing return Vector{Node}([]) end
    return val
end

function gather_hyb_descendants_outside_of_path_recur(node::Node, path::Array{Node}; stopathyb::Bool=false)
    children = getchildren(node)
    if length(children) == 0 && !(node in path)
        if !node.hybrid return nothing end  # nothing --> entire call will result in FALSE, so quit this loop quickly
        return Vector{Node}([node])
    end
    if stopathyb && node.hybrid return Vector{Node}([node]) end

    descendants = Vector{Node}([])
    for child in children
        if child in path continue end
        
        recur_descendants = gather_hyb_descendants_outside_of_path_recur(child, path, stopathyb=stopathyb)
        if recur_descendants === nothing
            return nothing
        end
        descendants = vcat(descendants, recur_descendants)
    end
    return descendants
end


TIEWARNING = false
"""
    findoptQ(D::Matrix{Float64}, idxpairs::Vector{Tuple{<:Integer, <:Integer}})

Finds the minimizer (i*, j*) among all pairs (i, j) in idxpairs for Q, a matrix computed from D.
"""
function findoptQidx(D::AbstractMatrix{Float64}, validpairs::BitArray, compat_trees::AbstractVector{HybridNetwork}; namelist=nothing, use_heuristic::Bool=true, max_sorted_entries::Int64=5)
    # max_sorted_entries: push!(sorted_values, ...) is the majority of algorithm run time. In an effort to combat this,
    #                     we iteratively cap the maximum size of sorted_values, b/c we usually don't actually need to
    #                     evaluate every `qij`, just the smallest 1 or 2
    global TIEWARNING

    n = size(D)[1]
    sums = sum(D, dims=1)
    sorted_values = SortedSet{Tuple{Float64, Tuple{Int64, Int64}}}()
    max_val = (-Inf, (-Inf, -Inf))

    # New, faster
    for i = 1:n
        for j = (i+1):n
            if !validpairs[i, j] continue end
            # TODO: refactor code so that we don't have to re-calculate `D`
            #       every time and specific entries are instead refactored when changed,
            #       then this function just receives the already sorted list of qij's
            #       instead of the matrix `D` itself

            qij = (n-2) * D[i,j] - sums[i] - sums[j]

            # push!(sorted_values, (qij, (i, j)))
            if length(sorted_values) < max_sorted_entries || qij < max_val[1]
                if length(sorted_values) >= max_sorted_entries
                    delete!(sorted_values, max_val)
                end

                push!(sorted_values, (qij, (i, j)))
                max_val = maximum(sorted_values)
            end
        end
    end

    if length(sorted_values) == 0
        if max_sorted_entries < n
            return findoptQidx(D, validpairs, compat_trees, max_sorted_entries=2*max_sorted_entries, namelist=namelist, use_heuristic=use_heuristic)
        end
        throw(SolutionDNEError())
    end

    if !use_heuristic
        return first(sorted_values)[2]
    else
        for (_, (i, j)) in sorted_values
            if are_compatible_after_merge(compat_trees, namelist[i], namelist[j])
                return (i, j)
            end
        end

        if max_sorted_entries < n
            return findoptQidx(D, validpairs, compat_trees, max_sorted_entries=2*max_sorted_entries, namelist=namelist, use_heuristic=use_heuristic)
        end
        throw(ErrorException("No compatible merge found."))
    end
end


"""
    findvalidpairs(constraints::Vector{HybridNetwork}, constraint_sibling_pairs, namelist::AbstractVector{<:AbstractString})

Finds all valid sibling pairs among the constraint networks.
"""
function findvalidpairs(constraints::Vector{HybridNetwork}, constraint_sibling_pairs, namelist::AbstractVector{<:AbstractString})
    n = length(namelist)

    # Shorthand functions for name lookups (we'll be doing a lot of these if there are many constraints)
    namedict = Dict{AbstractString, Int64}([name => i for (i, name) in enumerate(namelist)])
    idx(name::AbstractString) = namedict[name]

    # initialize matrix
    validpairs = falses(n, n, 2)    # [i, j, 1]: false if pair not seen together yet
                                    # [i, j, 2]: false = invalid pair, true = valid pair

    # go through the constraint networks and validate/invalidate pairs
    for (net_idx, net) in enumerate(constraints)
        if net.numtaxa == 1 continue end
        leafidxs = [idx(leaf.name) for leaf in net.leaf]

        # What is up with this xor (⊻)?
        # If we have seen a given entry before and NOT found that pair to be siblings (1 0) we want that to stay the same (1 stays 1)
        # if we have seen an entry AND FOUND that pair to be siblings (1 1) we want to change to (0 0) so it can be overridden
        # that's what this does.
        validpairs[leafidxs, leafidxs, 1] .= xor.(validpairs[leafidxs, leafidxs, 1], validpairs[leafidxs, leafidxs, 2])
        validpairs[leafidxs, leafidxs, 2] .= false

        # Find valid sibling pairs
        nodepairs = constraint_sibling_pairs[net_idx]
        # returned as a BitArray w/ indices mirroring net.leaf, need to convert to idxs
        nodestoidx(nodepair) = CartesianIndex(idx(nodepair[1].name), idx(nodepair[2].name))
        pairidxs = map(nodestoidx, nodepairs)

        # Set valid pair idxs to 1 only if they are either 1 or -1
        # if a pair idx is 0 then it was invalid elsewhere and needs to stay 0
        for idx in pairidxs
            if !validpairs[idx[1], idx[2], 1]
                # Enter twice, once in the upper triangle and once in the lower+
                validpairs[idx[1], idx[2], 2] = validpairs[idx[2], idx[1], 2] = true
            end
        end

        # Mark all these indices as seen together now
        validpairs[leafidxs, leafidxs, 1] = fill!(validpairs[leafidxs, leafidxs, 1], true)
    end

    # Any pairs that still have not been seen are valid
    return validpairs[:, :, 2] .|| .!validpairs[:, :, 1]
end
findvalidpairs(net::HybridNetwork, namelist::AbstractVector{<:AbstractString}; kwargs...) = findvalidpairs([net], namelist; kwargs...)


"""
    findsiblingpairs(net::HybridNetwork)

Finds sibling pairings for all the leaves in a single network.
These pairs are valid for `net` but may not be valid when the
    other constraint networks are also considered.
Returns a vector of tuples of nodes corresponding to siblings.
"""
function findsiblingpairs(net::HybridNetwork; force_unrooted::Bool=true, kwargs...)
    pairs = BitArray(undef, net.numtaxa, net.numtaxa)
    pairs .= false

    # new, simpler method just using Graphs.jl
    # hybedges = Vector{Any}([edge for edge in net.edge if (edge.hybrid && !edge.ismajor)])
    # push!(hybedges, nothing)
    
    node_map = Dict(leaf => i for (i, leaf) in enumerate(net.leaf))

    for nodei_idx = 1:net.numtaxa
        nodei = net.leaf[nodei_idx]
        if nodei == getroot(net) continue end
        neighbors = [child for child in getchildren(getparent(nodei)) if child.leaf && !(child == nodei)]

        if getparent(nodei) == getroot(net) && force_unrooted
            for child in getchildren(getparent(nodei))
                for child_child in getchildren(child)
                    if child_child.leaf && !(child_child == nodei)
                        push!(neighbors, child_child)
                    end
                end
            end
        end

        for neighbor in neighbors
            idx1 = min(nodei_idx, node_map[neighbor])
            idx2 = max(nodei_idx, node_map[neighbor])
            pairs[idx1, idx2] = true
        end
    end
    
    # Convert pairs to Tuples of Nodes
    nodepairs = Array{Tuple{Node, Node}}(undef, sum(pairs))
    for (arr_idx, leaf_idxs) in enumerate(findall(pairs))
        nodepairs[arr_idx] = (net.leaf[leaf_idxs[1]], net.leaf[leaf_idxs[2]])
    end

    return nodepairs
end


"""
Helper function - returns true if the node `node` is connected directly to
the root of `net`, or if the path connecting `node` to `net`'s root is
multiple redundant edges (e.g. returns true for A in the net: "((((((A))))), (B,C));")
"""
function has_direct_root_connection(net::HybridNetwork, node::Node)
    curr = node
    last_node = nothing

    while true
        if curr == getroot(net) break end
        if length(getparents(curr)) > 1 return false end
        
        children = getchildren(curr)
        for child in children
            if child != curr && child.hybrid
                try
                    if curr == getparent(getparentedge(child)) return false end
                catch
                end
            end
        end

        if length(getparents(curr)) == 0 throw(ErrorException("Unexpected case.")) end

        last_node = curr
        curr = getparents(curr)[1]
    end
    return true
end


"""

Helper function for manipulating graph adjacency lists.
"""
function vecreplace!(vec::Vector{T}, pattern::T, with::T) where T
    for (i, val) in enumerate(vec)
        if val == pattern
            vec[i] = with
        end
    end
end



getnodes(n::Node) = reduce(vcat, [[child for child in e.node if child != n] for e in n.edge], init=Vector{Node}())


"""
Helper function
"""
function log_edge_path_retics_from_node(
    node::Node, edgesinpath::AbstractVector{Edge}, relevanttoi::Bool, hybedgelogged::Bool,
    hybedge::Union{Nothing, Edge}, reticmap::ReticMap, subnetedgei::Edge, subnetedgej::Edge,
    net::HybridNetwork, edgeinpath_logged::Bool, nodesinpath::AbstractVector{Node},
    nodei::Node, nodej::Node
)
    node_edges = node.edge
    isminorhyb = [e.hybrid && !e.ismajor for e in node_edges]

    if sum(isminorhyb) == 1
        edge = node_edges[isminorhyb][1]

        if !(edge in edgesinpath)
            # Edge is not traversed in the path
            fromorto = ifelse(getchild(edge) == node, "to", "from")
            
            if !hybedgelogged || edge != hybedge
                trylogretic_single!(reticmap, edge, ifelse(relevanttoi, subnetedgei, subnetedgej), fromorto)
                if edge == hybedge hybedgelogged = true end

                # Special edge case: constraint only has 2 leaves (i.e. this is its last merge)
                #                    and the minor edge is OUTSIDE the path, then we need to log
                #                    the other half of it now
                if length(net.leaf) == 2
                    # We use trylogretic_single! b/c this retic may have already been logged previously in this iteration
                    trylogretic_single!(reticmap, edge, ifelse(relevanttoi, subnetedgej, subnetedgei), ifelse(fromorto == "from", "to", "from"))
                end
            end
        else
            # Edge *IS* traversed in the path
            if !edgeinpath_logged
                if getchild(edge) == node
                    # edge goes from j --> i
                    logretic!(reticmap, edge, subnetedgej, "from")
                    logretic!(reticmap, edge, subnetedgei, "to")
                else
                    # edge goes from i --> j
                    logretic!(reticmap, edge, subnetedgei, "from")
                    logretic!(reticmap, edge, subnetedgej, "to")
                end

                edgeinpath_logged = true
            end
        end

    elseif sum(isminorhyb) == 2
        hyb_edges = node_edges[isminorhyb]
        if hyb_edges[1] in edgesinpath && hyb_edges[2] in edgesinpath
            error("This scenario should not be possible. Please submit an issue on GitHub.")
        end

        # `edge_i` corresponds to `node_i`, and `_j` to `_j`.
        # only one of these edges will appear in `edgesinpath`.
        # say `hyb_edges[1]` (call it `hybedge1`) appears in `edgesinpath`.
        # if its index in `edgesinpath` is the index of `node` in `nodesinpath`
        # MINUS ONE, then the hybrid is directed from `nodej` into `nodei`.
        # Otherwise, it is directed from `nodei` into `nodej`.
        inpath_edge = ifelse(hyb_edges[1] in edgesinpath, hyb_edges[1], hyb_edges[2])
        notinpath_edge = ifelse(hyb_edges[1] in edgesinpath, hyb_edges[2], hyb_edges[1])

        inpath_edge_childnode = getchild(inpath_edge)
        inpath_edge_parentnode = getparent(inpath_edge)
        inpath_from_i_to_j = findfirst([inpath_edge_parentnode] .== nodesinpath) < findfirst([inpath_edge_childnode] .== nodesinpath)

        notinpath_direction, notinpath_subnet = find_not_in_path_edge_vector(notinpath_edge, subnetedgei, subnetedgej, nodesinpath, nodei, nodej)

        # Now we know the direction of the retic, so let's place them.
        if inpath_from_i_to_j
            logretic!(reticmap, inpath_edge, subnetedgei, "from")
            trylogretic!(reticmap, inpath_edge, subnetedgej, "to")
            logretic!(reticmap, notinpath_edge, notinpath_subnet, notinpath_direction)
        else
            logretic!(reticmap, inpath_edge, subnetedgej, "from")
            trylogretic!(reticmap, inpath_edge, subnetedgei, "to")
            logretic!(reticmap, notinpath_edge, notinpath_subnet, notinpath_direction)
        end
        edgeinpath_logged = true
    end


    # Gather all of the descendant tips of the current node. If all of these nodes *except*
    # those in `nodesinpath` are hybrids *that don't lead anywhere* (i.e. their "to" portion
    # has been logged), then log the "from" portion of each of those retics and remove the edges

    if length(net.hybrid) > 0
        node_descendants = gather_hyb_descendants_outside_of_path(node, nodesinpath)
        if length(node_descendants) > 0 && all(d.hybrid && length(getchildren(d)) == 0 for d in node_descendants)
            relevant_subnet_edge = relevanttoi ? subnetedgei : subnetedgej
            for hyb in node_descendants
                trylogretic_single!(reticmap, hyb, relevant_subnet_edge, "from")
                # -- this might need to be changed to `trylogretic!`, not sure --
            end
        end

        hyb_descendants = gather_hyb_descendants_outside_of_path(node, nodesinpath, stopathyb=true)
        if length(hyb_descendants) > 0 && all(d.hybrid for d in hyb_descendants)
            relevant_subnet_edge = relevanttoi ? subnetedgei : subnetedgej
            for hyb in hyb_descendants
                trylogretic_single!(reticmap, hyb, relevant_subnet_edge, "from")
            end
        end
    end

    return edgeinpath_logged
end


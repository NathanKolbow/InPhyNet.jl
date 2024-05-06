
ERROR_DIR = "/mnt/dv/wid/projects4/SolisLemus-network-merging/error_logs/"
function load_next_debug_data()
    global ERROR_DIR
    if length(readdir(ERROR_DIR)) == 0
        error("No more error logs.")
    end

    error_id = readdir(ERROR_DIR)[1]
    error_id = split(error_id, "dmat_")[2]
    error_id = split(error_id, ".csv")[1]

    noisy_D = Matrix(CSV.read(joinpath(ERROR_DIR, "dmat_$(error_id).csv"), DataFrame))

    log_file = "error_$(error_id).log"
    log_lines = readlines(joinpath(ERROR_DIR, log_file))
    
    # Grab the constraints
    start_idx = 3
    end_idx = 3
    while log_lines[end_idx] != ""
        end_idx += 1
        if end_idx > length(log_lines)
            error("$(log_file) has unknown format.")
        end
    end
    end_idx -= 1

    constraints = [readTopology(c) for c in log_lines[start_idx:end_idx]]
    truenet = readTopology(log_lines[end_idx+3])
    true_D, namelist = majorinternodecount(truenet)

    return truenet, constraints, noisy_D, true_D, namelist, parse(Int64, error_id)
end


function step_inphynet_starter_vars(D, constraints, namelist)
    D_iter = deepcopy(D)
    cs_iter = deepcopy(constraints)
    namelist_iter = deepcopy(namelist)
    subnets = Vector{SubNet}([SubNet(i, namelist[i]) for i in 1:length(namelist)])
    reticmap = reticmap = ReticMap(cs_iter)
    rootretics = Array{Union{Tuple, Edge, Nothing}}(undef, length(constraints))
    for (i, c) in enumerate(constraints)
        hybridbools = [edge.hybrid for edge in c.node[c.root].edge]
        if any(hybridbools)
            sum(hybridbools) == 1 || error("Two reticulations coming out of root, this is not accounted for!")
            hybedge = c.node[c.root].edge[hybridbools][1]
            rootretics[i] = hybedge
        else
            # None of the nodes coming out of the root are hybrids, but now
            # we need to make sure that none of the nodes coming out of the
            # root lead to _only_ hybrids, thus effectively making the root
            # lead directly to hybrids.
            children = getchildren(c.node[c.root])
            hybridbools = [[edge.hybrid for edge in child.edge] for child in children]

            @show hybridbools
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

    return D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed
end


function step_inphynet!(D, constraints, namelist, subnets, reticmap, rootretics, rootreticprocessed)
    major_tree_only = false

    n = size(D, 1)

    ## SINGLE ALGO STEP ##
    possible_siblings = InPhyNet.findvalidpairs(constraints, namelist, major_tree_only = major_tree_only)
    
    # Find optimal (i, j) idx pair for matrix Q
    i, j = InPhyNet.findoptQidx(D, possible_siblings)
    @show (namelist[i], namelist[j])

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
    ######################
    return D, constraints, namelist, subnets, reticmap, rootretics, rootreticprocessed
end


function find_problematic_constraints(D, constraints, namelist)
    prob_constraints = Vector{HybridNetwork}([])
    
    for c in constraints
        try
            netnj(D, [c], namelist)
        catch e
            if !(typeof(e) <: SolutionDNEError)
                push!(prob_constraints, c)
            end
        end
    end
    return prob_constraints
end


function reduce_D_namelist(D, constraints, namelist)
    keep_names = reduce(vcat, [[leaf.name for leaf in c.leaf] for c in constraints])
    remove_names = setdiff(namelist, keep_names)
    for name in remove_names
        idx = findfirst(namelist .== name)
        idx_filter = setdiff(1:length(namelist), idx)
        D = D[idx_filter, idx_filter]
        namelist = namelist[idx_filter]
    end
    return D, namelist
end


function remove_error_files(error_id)
    rm(joinpath(ERROR_DIR, "dmat_$(error_id).csv"))
    rm(joinpath(ERROR_DIR, "error_$(error_id).log"))
end
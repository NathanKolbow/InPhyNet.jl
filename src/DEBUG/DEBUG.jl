using Revise
using InPhyNet, PhyloNetworks, CSV, DataFrames

include("/mnt/dv/wid/projects4/SolisLemus-network-merging/src/DEBUG/debug-helpers.jl")


# Get rid of non-reproducible errors
while true
    truenet, constraints, D, true_D, namelist, error_id = load_next_debug_data()
    try
        netnj(D, constraints, namelist)
        @info "Merge successful in $(error_id), removing."
        remove_error_files(error_id)
    catch e
        if typeof(e) <: SolutionDNEError
            @info "SolutionDNEError in $(error_id), removing."
            remove_error_files(error_id)
        else
            @info "Found case for netid n$(truenet.numTaxa-1)r$(truenet.numHybrids) ($(typeof(e)), $(error_id))"
            break
        end
    end
end

# Load error data
truenet, constraints, D, true_D, namelist, error_id = load_next_debug_data()
@info "Error ID: $(error_id)"

# Confirm that an error occurs
netnj(D, constraints, namelist)

# Narrow down which constraint(s) cause the error
cs = find_problematic_constraints(D, constraints, namelist)
cs = constraints[3:4]
D_reduced, namelist_reduced = reduce_D_namelist(D, cs, namelist)    # sometimes the problem will still occur when D & namelist are reduced, sometimes not

#### DEBUG ####

D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed = step_inphynet_starter_vars(D, cs, namelist)
while size(D_iter, 1) > 1
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
        step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
end

mnet = InPhyNet.HybridNetwork(subnets[1].nodes, subnets[1].edges)
mnet.root = mnet.numNodes
mnet.node[mnet.root].name = "root"
mnet = InPhyNet.placeretics!(mnet, reticmap, copy_retic_names=false)
InPhyNet.removeplaceholdernames!(mnet)




# MWE
net = readTopology("(((A)#H1,B),#H1);")
r = InPhyNet.ReticMap([net])

nodei = net.leaf[1]
nodej = net.leaf[2]

graph = InPhyNet.Graph(net, includeminoredges=true)

idxnodei = findfirst(net.node .== [nodei])
idxnodej = findfirst(net.node .== [nodej])
edgepath = a_star(graph, idxnodei, idxnodej)

nodesinpath = Array{Node}(undef, length(edgepath)+1)
edgesinpath = Array{Edge}(undef, length(edgepath))
for (i, gedge) in enumerate(edgepath)
    srcnode = net.node[gedge.src]
    dstnode = net.node[gedge.dst]

    if i == 1 nodesinpath[1] = net.node[gedge.src] end
    nodesinpath[i+1] = dstnode

    netedge = filter(e -> (srcnode in e.node) && (dstnode in e.node), dstnode.edge)[1]
    edgesinpath[i] = netedge
end

# find the node that should be the new tip after merging
newtip = nodesinpath[1]
while length(getparents(newtip)) > 0 && getparent(newtip) in nodesinpath
    newtip = getparent(newtip)
end





InPhyNet.updateconstraints!("A", "B", [net], r, Edge(-1), Edge(-1))
sum(r.map[collect(keys(r.map))[1]] .!== nothing) > 1
###############






# Remove the error files once resolved
remove_error_files(error_id)






using PhyloPlots
PhyloPlots.plot(constraints[1])





cc = readTopology("((A,(B,#H1)),(C)#H1);")
_, namelist = majorinternodecount(cc)
D = Matrix{Float64}([0. 2. 3.; 2. 0. 0.; 3. 0. 0.])
mnet = netnj(D, [cc], namelist)
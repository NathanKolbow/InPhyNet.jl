ENV["JULIA_DEBUG"] = Main
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
cs = [constraints[14]]
D_reduced, namelist_reduced = reduce_D_namelist(D, cs, namelist)    # sometimes the problem will still occur when D & namelist are reduced, sometimes not
netnj(D_reduced, cs, namelist_reduced)


for c in cs for node in c.node
    if node.name == "" node.name = "__$(node.number)__" end
end end

#### DEBUG ####
D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
    # step_inphynet_starter_vars(D_reduced, cs, namelist_reduced)
    step_inphynet_starter_vars(D, cs, namelist)
    # step_inphynet_starter_vars(D, constraints, namelist)

i = 1

ff() = findfirst(node.name == "H39" for node in cs_iter[1].node)

non_nothing_vals(rmap) = sum(sum(rmap.map[key] .!== nothing) for key in keys(rmap.map))
while true
    @show i
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
        step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
    i += 1
end

# (391, 378)
for i=1:1
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
        step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
end

graph, W, nodesinpath, edgesinpath = InPhyNet.find_valid_node_path(cs_iter[1], "t383", "t377")


# (t403, t422)
idxnodei = findfirst(node.name == "t403" for node in cs_iter[1].node)
idxnodej = findfirst(node.name == "t422" for node in cs_iter[1].node)
nodei = cs_iter[1].node[idxnodei]
nodej = cs_iter[1].node[idxnodej]

graph, W = InPhyNet.Graph(cs_iter[1], includeminoredges=true, withweights=true, minoredgeweight=1.51, removeredundantedgecost=true)
edgepath = a_star(graph, idxnodei, idxnodej, W)

nodesinpath = Array{Node}(undef, length(edgepath)+1)
edgesinpath = Array{Edge}(undef, length(edgepath))
for (i, gedge) in enumerate(edgepath)
    srcnode = cs_iter[1].node[gedge.src]
    dstnode = cs_iter[1].node[gedge.dst]

    if i == 1 nodesinpath[1] = cs_iter[1].node[gedge.src] end
    nodesinpath[i+1] = dstnode

    netedge = filter(e -> (srcnode in e.node) && (dstnode in e.node), dstnode.edge)[1]
    edgesinpath[i] = netedge
end




# (t305, t301)
for i=1:1
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
        step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
end

non_nothing_vals(rmap) = sum(sum(rmap.map[key] .!== nothing) for key in keys(rmap.map))
while non_nothing_vals(reticmap) < 1
    @show i
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
        step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
    i += 1
end

while i < 69
    @show i
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
        step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
    i += 1
end

for node in cs_iter[1].node
    if node.name == "" node.name = "__$(node.number)__" end
end

graph, W = InPhyNet.Graph(cs_iter[1], includeminoredges=true, withweights=true, minoredgeweight=1.51, removeredundantedgecost=true)
nodei = cs_iter[1].node[findfirst([node.name == "t66" for node in cs_iter[1].node])]
nodej = cs_iter[1].node[findfirst([node.name == "t63" for node in cs_iter[1].node])]
idxnodei = findfirst(cs_iter[1].node .== [nodei])
idxnodej = findfirst(cs_iter[1].node .== [nodej])
edgepath = a_star(graph, idxnodei, idxnodej, W)

# (t66, t63)
for i=1:1
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
        step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
end





# MWE
net = readTopology(writeTopology(temp[1]))
# net = readTopology(writeTopology(cs[1]))

D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed = step_inphynet_starter_vars(D_reduced, [net], namelist_reduced)

for i=1:4
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
    step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
end

#
D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
    step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
###############






# Remove the error files once resolved
remove_error_files(error_id)






using PhyloPlots
PhyloPlots.plot(constraints[1])





cc = readTopology("((A,(B,#H1)),(C)#H1);")
_, namelist = majorinternodecount(cc)
D = Matrix{Float64}([0. 2. 3.; 2. 0. 0.; 3. 0. 0.])
mnet = netnj(D, [cc], namelist)
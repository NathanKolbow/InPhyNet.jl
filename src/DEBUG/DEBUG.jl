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
D_reduced, namelist_reduced = reduce_D_namelist(D, cs, namelist)    # sometimes the problem will still occur when D & namelist are reduced, sometimes not
netnj(D_reduced, cs, namelist_reduced)


temp = [pruneTruthFromDecomp(cs[1], ["t6", "t12", "t11", "t10", "t4", "t2", "t1", "t3", "t5"])]
D_reduced, namelist_reduced = reduce_D_namelist(D, temp, namelist)
netnj(D_reduced, temp, namelist_reduced)

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
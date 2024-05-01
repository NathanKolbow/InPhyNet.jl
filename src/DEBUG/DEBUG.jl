using Revise
using InPhyNet, PhyloNetworks, CSV, DataFrames

include("debug-helpers.jl")


# Get rid of non-reproducible errors
while true
    truenet, constraints, D, true_D, namelist, error_id = load_next_debug_data()
    try
        netnj(D, constraints, namelist)
        @info "Merge successful in $(error_id), removing."
        remove_error_files(error_id)
    catch e
        if typeof(e) <: SolutionDNEError
            println("SolutionDNEError in $(error_id), removing.")
            remove_error_files(error_id)
        else
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

#### DEBUG ####

D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed = step_inphynet_starter_vars(D, cs, namelist)
while cs_iter[1].numTaxa > 1
    D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed =
        step_inphynet!(D_iter, cs_iter, namelist_iter, subnets, reticmap, rootretics, rootreticprocessed)
end

mnet = InPhyNet.HybridNetwork(subnets[1].nodes, subnets[1].edges)
mnet.root = mnet.numNodes
mnet.node[mnet.root].name = "root"
mnet = InPhyNet.placeretics!(mnet, reticmap, copy_retic_names=false)
InPhyNet.removeplaceholdernames!(mnet)


###############






# Remove the error files once resolved
remove_error_files(error_id)






using PhyloPlots
PhyloPlots.plot(constraints[1])





cc = readTopology("((A,(B,#H1)),(C)#H1);")
_, namelist = majorinternodecount(cc)
D = Matrix{Float64}([0. 2. 3.; 2. 0. 0.; 3. 0. 0.])
mnet = netnj(D, [cc], namelist)
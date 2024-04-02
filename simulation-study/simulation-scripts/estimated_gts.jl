# Full pipeline going from truenet > simulated gene trees > simulated sequences > inferred gene trees > inferred networks > netmerge

if length(ARGS) != 5
    error("Usage: julia estimated_gts.jl \"<true network abbreviation>\" <replicate number> <maximum subset size> \"<distance method>\", <number of gene trees>")
end

###### Input parsing ######
netid = ARGS[1]
replicatenum = parse(Int64, ARGS[2])
maxsubsetsize = parse(Int64, ARGS[3])
dmethod = ARGS[4]
ngt = parse(Int64, ARGS[5])
###########################

include("helpers/helpers.jl")

# 1. gather ground truth network, constraint, distance matrix, and namelist
println("Loading data")
truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, dmethod)

# 1.1 set seed
seed = parse(Int64, "$(truenet.numTaxa)42$(truenet.numHybrids)42$(replicatenum)")

# 2. simulate gene trees
println("Simulating gene trees")
using PhyloCoalSimulations
truegt_file = "./data/truegt_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).treefile"
gts = nothing

if isfile(truegt_file)
    gts = readMultiTopology(truegt_file)
else
    Random.seed!(seed)
    gts = simulatecoalescent(truenet, ngt, 1)
    writeMultiTopology(gts, truegt_file)
end

# 3. simulate sequences w/ seq-gen
println("Simulating sequences")
seq_file = "./data/seqfile_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).fasta"
if !isfile(seq_file)
    Random.seed!(seed)
    # simulate sequences here
end

# 4. infer gene trees w/ iqtree
println("Inferring gene trees")
estgt_file = "./data/estgt_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).treefile"
if !isfile(estgt_file)
    Random.seed!(seed)
    # infer gene trees here
end

# 5. infer networks w/ SNaQ
println("Inferring networks")
net_file = "./data/estnets_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).netfile"
if !isfile(estgt_file)
    Random.seed!(seed)
    # infer networks here
end

# 6. InPhyNet inference
println("Merging networks")
est_constraints = readMultiTopology(net_file)
est_gts = readMultiTopology(estgt_file)
est_D, est_namelist = calculateAGIC(est_gts)

mnet = netnj(est_D, est_constraints, est_namelist)
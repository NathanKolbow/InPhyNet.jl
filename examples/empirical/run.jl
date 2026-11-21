# Empirical example: InPhyNet on the Lepidoptera phylogenomic data of
# Kawahara et al. (2019, PNAS, doi:10.1073/pnas.1907847116).
# Reproduces the results shown on the "Empirical example" docs page.
# This script is NOT run during the docs build.
#
# Data: download Supplementary Archive 5 from Dryad (doi:10.5061/dryad.j477b40),
# unpack it, and set `datadir` below to the folder that contains
# "Supplementary Archive 5". The gene trees are in
# Supplementary Archive 5/SA5_aa_ml_trees/ (1,930 files ending in .treefile).
#
# Run from the repo root, e.g. with 10 worker processes:
#   julia --project=examples/empirical -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=examples/empirical -p 10 examples/empirical/run.jl
# The Julia and package versions used are written to `outdir`/versions.txt.

const datadir = "/path/to/data"                           # <-- set this
const outdir  = joinpath(datadir, "inphynet-empirical")   # results go here
const minsize, maxsize = 5, 15   # subset sizes for centroid edge decomposition
const hmax  = 1                  # reticulations per constraint network
const nruns = 10                 # SNaQ runs per subset
const seed  = 2026
const outgroup = "Coleoptera_Tenebrionidae_Tenebrioninae_Tribolium_castaneum_TCAST"

using Distributed
@everywhere using PhyloNetworks, SNaQ
using InPhyNet, DataFrames, Dates, Pkg, PhyloPlots, RCall
mkpath(joinpath(outdir, "snaq"))
timings = Dict{String,Float64}()

# 1. Read the gene trees
gtdir = joinpath(datadir, "Supplementary Archive 5", "SA5_aa_ml_trees")
gtfiles = filter(endswith(".treefile"), readdir(gtdir))
gts = [readnewick(joinpath(gtdir, f)) for f in gtfiles]
@info "Read $(length(gts)) gene trees"

# 2. Distance matrix: average gene-tree internode distance
timings["AGID"] = @elapsed D, namelist = calculateAGID(gts)
@info "$(length(namelist)) taxa"

# 3. Guide tree (neighbor joining) and subsets
guide = nj(DataFrame(D, namelist))
writenewick(guide, joinpath(outdir, "guide_tree.tre"))
subsets = centroid_edge_decomposition(guide, minsize, maxsize)
open(joinpath(outdir, "subsets.txt"), "w") do io
    for s in subsets; println(io, join(sort(s), ",")); end
end
@info "$(length(subsets)) subsets, sizes $(sort(length.(subsets)))"

# 4. Constraint networks: SNaQ on each subset (its runs are spread over the workers)
constraints = Vector{HybridNetwork}(undef, length(subsets))
for (i, subset) in enumerate(subsets)
    prefix = joinpath(outdir, "snaq", "subset$(i)")
    if isfile(prefix * ".out")   # resume: this subset already finished
        constraints[i] = readsnaqnetwork(prefix * ".out")
        continue
    end
    @info "SNaQ on subset $i/$(length(subsets)) ($(length(subset)) taxa)"
    t = @elapsed begin
        # keep gene trees with at least 4 taxa of the subset (fewer give no quartets)
        sub_gts = filter(g -> count(in(subset), tiplabels(g)) >= 4, gts)
        df = readtrees2CF(prune_networks(sub_gts, subset), writeTab=false, writeSummary=false)
        # starting tree: the guide tree pruned to the subset
        start = readnewick(writenewick(prune_network(guide, subset)))
        constraints[i] = snaq!(start, df, hmax=hmax, runs=nruns, seed=seed+i, filename=prefix)
    end
    timings["SNaQ subset $i ($(length(subset)) taxa)"] = t
end

# 5. Merge with InPhyNet and root on the outgroup
timings["InPhyNet merge"] = @elapsed net = inphynet(D, constraints, namelist)
rootatnode!(net, outgroup)
writenewick(net, joinpath(outdir, "inphynet_network.tre"))
@info "Merged network: $(net.numtaxa) taxa, $(net.numhybrids) reticulations"

# 6. Plot, with tip labels shortened to "Genus species (Superfamily)"
function shortlabel(s)
    p = split(s, "_")   # "Superfamily_Family_[Subfamily_]Genus_species_ID"
    for i in 2:length(p)-1
        if occursin(r"^[A-Z][a-z]+$", p[i]) && !occursin(r"(idae|inae|oidea)$", p[i]) &&
           occursin(r"^[a-z]+$", p[i+1])
            return string(p[i], " ", p[i+1], " (", p[1], ")")
        end
    end
    return string(p[1], ": ", join(p[2:min(3, end)], " "))  # no species name
end
plotnet = deepcopy(net)
for l in plotnet.leaf; l.name = shortlabel(l.name); end
R"pdf(NULL)"; r0 = plot(plotnet; tipcex=0.55, showgamma=false); R"dev.off()"
R"png"(joinpath(outdir, "empirical_network.png"), width=1600, height=5200, res=150)
R"par"(mar=[0,0,0,0])
plot(plotnet; tipcex=0.55, showgamma=false, xlim=(r0.xmin, r0.xmin + 1.3*(r0.xmax - r0.xmin)))
R"dev.off()"

# 7. Record versions and timings
open(joinpath(outdir, "versions.txt"), "w") do io
    println(io, "date: ", now())
    println(io, "julia: ", VERSION)
    println(io, "workers: ", nworkers())
    for (_, p) in Pkg.dependencies()
        p.is_direct_dep && println(io, p.name, ": ", p.version)
    end
    println(io, "\ntimings (seconds):")
    for (k, v) in sort(collect(timings)); println(io, "  ", k, ": ", round(v, digits=1)); end
end

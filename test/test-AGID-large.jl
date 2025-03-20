using StatsBase, Random, PhyloNetworks, Test
using InPhyNet

est_gts = readmultinewick(joinpath(@__DIR__, "examples/Best.33.FAA.tre"))[21:50];
D, namelist = calculateAGID(est_gts);

Random.seed!(42)
test_pairs = [sample(namelist, 2, replace=false) for j=1:10]
manual_dists = zeros(10) .- 1

for (j, (t1, t2)) in enumerate(test_pairs)
    _i1 = findfirst(t_name -> t_name == t1, namelist)
    _i2 = findfirst(t_name -> t_name == t2, namelist)
    D_val = D[_i1, _i2]

    manual_dists[j] = 0.0
    n = 0
    for gt in est_gts
        if !(t1 in tiplabels(gt)) || !(t2 in tiplabels(gt)) continue end
        n += 1

        iter_D, iter_namelist = internodedistance(gt)
        _i1 = findfirst(t_name -> t_name == t1, iter_namelist)
        _i2 = findfirst(t_name -> t_name == t2, iter_namelist)
        manual_dists[j] += iter_D[_i1, _i2]
    end
    manual_dists[j] /= n

    @test manual_dists[j] == D_val
end

using StatsBase, Random

est_gts = readMultiTopology("tests/examples/Best.33.FAA.tre");
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
        if !(t1 in tipLabels(gt)) || !(t2 in tipLabels(gt)) continue end
        n += 1

        iter_D, iter_namelist = internodedistance(gt)
        _i1 = findfirst(t_name -> t_name == t1, iter_namelist)
        _i2 = findfirst(t_name -> t_name == t2, iter_namelist)
        manual_dists[j] += iter_D[_i1, _i2]
    end
    manual_dists[j] /= n

    @test manual_dists[j] == D_val

    # @info "$(D_val) =?= $(manual_dists[j]): $(D_val == manual_dists[j])"

end
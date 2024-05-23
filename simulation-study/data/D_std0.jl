using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")

out_path = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/data/D_std0.csv"
out_df = DataFrame(
    net_id=String[],
    rep=Int64[],
    std0=Float64[]
)

# for net_id in get_all_net_ids()
for net_id in ["n50r2", "n100r10"]
    @info net_id
    std0s = zeros(100) .- 1.

    Threads.@threads for rep=1:100
        # truenet, _, _, _ = loadPerfectData(net_id, rep, 10, "AGIC")
        truenet = readMultiTopology(getNetworkFilepath(net_id))[rep]
        avg_bl = get_avg_bl(truenet)
        newick = writeTopology(truenet)
        newick = "($(newick[1:(length(newick)-1)]):$(avg_bl),OUTGROUP:1.0);"
        truenet = readTopology(newick)

        D, _ = majorinternodecount(truenet)
        std0 = upperTriangStd(D)
        std0s[rep] = std0
    end

    for i=1:100
        push!(out_df, [net_id, i, std0s[i]])
    end
end

CSV.write(out_path, out_df)

include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")

function test_retry(net_id::String, rep_num::Int64, m::Int64, max_retries::Int64, nsim::Int64=100)
    # 1. Gather true data
    truenet, constraints, D, namelist = loadPerfectData(net_id, rep_num, m, "AGIC")

    no_retry_successes = AtomicCounter(0)
    retry_successes = AtomicCounter(0)

    Threads.@threads for i=1:nsim
        # 2. Generate noise
        std0 = upperTriangStd(D)
        gaussSd = rand(Uniform(0.25*std0, 1.5*std0))

        totalnnimoves = Int64(round(rand(Uniform(0, length(constraints) * 4))))
        nnimoves = sample(1:length(constraints), totalnnimoves, replace=true)
        nnimoves = Vector{Int64}([sum(nnimoves .== i) for i=1:length(constraints)])

        # 3. Copy onto data
        D_iter = deepcopy(D)
        addnoise!(D_iter, Normal(gaussSd, gaussSd))
        constraints_iter = copyConstraints(constraints)
        for (i, (c, nmoves)) in enumerate(zip(constraints_iter, nnimoves))
            for _=1:nmoves doRandomNNI!(c) end
            # constraintdiffs[i] = hardwiredClusterDistance(constraints[i], c)
        end

        # 4. Run netnj w/o retrying
        try
            mnet = netnj(D_iter, constraints_iter, namelist, supressunsampledwarning=true)
            @atomic :sequentially_consistent no_retry_successes.iterspassed += 1
        catch
        end

        # 5. Run netnj w/ retrying
        try
            mnet = silently() do
                InPhyNet.netnj_retry_driver(D_iter, constraints_iter, namelist, verbose=false, max_retry_attempts=max_retries)
            end
            if mnet !== nothing
                @atomic :sequentially_consistent retry_successes.iterspassed += 1
            end
        catch e
            rethrow(e)
            println("retry error")
        end
    end

    return no_retry_successes.iterspassed / nsim, retry_successes.iterspassed / nsim
end


results = []
for net_id in ["n50r2", "n100r10", "n200r10"]
    for m in [15, 20]
        noretry_rate, retry_rate = test_retry(net_id, 1, m, 100, 100)
        @info "$(net_id) #1, m=$(m) - $(noretry_rate) vs. $(retry_rate)"
        push!(restults, (net_id, m, noretry_rate, retry_rate))
    end
end

using PhyloCoalSimulations

function log(tag::String, msg::String, msg_color::Symbol=:cyan, tag_color::Symbol=:black)
    printstyled("[$(tag)] ", color=tag_color)
    printstyled("$(msg)\n", color=msg_color)
end


function print_basic_info()
    println("-------------------------")
    printstyled("$(nprocs()) CPUs\n", color=:red)
    printstyled("$(Threads.nthreads()) Threads\n", color=:green)
    printstyled("UTC $(Int64(round(1e9 * time())))\n", color=:yellow)
    println("-------------------------\n")
end


function simulate_gene_trees(truegt_file::String, ngt::Int64, seed)
    Random.seed!(seed)
    
    gts = nothing
    if isfile(truegt_file)
        gts = readMultiTopology(truegt_file)
        log("TRUE GT", "Loaded from file.")
    else
        Random.seed!(seed)
        log("TRUE GT", "Simulating true gene trees.")
        gts = simulatecoalescent(true_net, ngt, 1)
        writeMultiTopology(gts, truegt_file)
    end

    return gts
end


function est_gt_from_true_gt(gt::HybridNetwork, seq_file::String, estgt_file::String, data_dir::String, idx::Int64)
    # 1. simulate sequence w/ seq-gen (checkpointed)
    if !isfile(seq_file)
        log("SEQ", "Simulating sequences for #$(idx).")
        temp_gt_file::String = joinpath(data_dir, "temp_gt_$(abs(rand(Int64)))")
        writeTopology(gt, temp_gt_file)
        run_seqgen(0.036, temp_gt_file, seq_file)
        rm(temp_gt_file)
    else
        log("SEQ", "Sequences already simulated for #$(idx).")
    end

    # 2. infer gene trees w/ iqtree (checkpointed)
    if !isfile(estgt_file)
        log("IQ-TREE", "Running IQ-Tree on #$(idx).")
        run_iqtree(seq_file)
        mv("$(seq_file).treefile", estgt_file)

        clean_iqtree_files(seq_file)
        rm(seq_file)
    else
        log("IQ-TREE", "Results already exist for #$(idx).")
    end
end


function estimate_nj_tree(est_gts::AbstractVector{HybridNetwork})
    log("NJ", "Calculating AGIC.")
    est_D, est_namelist = calculateAGIC(est_gts)

    log("NJ", "Estimating NJ tree.")
    nj_df = DataFrame(est_D, est_namelist)
    nj_tre = nj(nj_df)

    return est_D, est_namelist, nj_tre
end


function subset_decomp(nj_tre::HybridNetwork, m::Int64)
    log("SUBSETS", "Decomposing into subsets.")
    subsets = sateIdecomp(nj_tre, m)

    if minimum([length(s) for s in subsets]) < 5
        log("SUBSETS", "Smallest subset must have at least 5 taxa for SNaQ.", :red, :red)
        log("SUBSETS", "Smallest subset size: $(minimum([length(s) for s in subsets]))", :red, :red)
        error("Smallest subset must have at least 5 taxa for SNaQ.")
        exit(100)
    end
    return subsets
end


function snaq_constraints(est_gts::AbstractVector{HybridNetwork}, output_file_prefix::String, subsets, true_net::HybridNetwork, nj_tre::HybridNetwork, seed::Int64)
    Random.seed!(seed)

    q, t = silently() do
        countquartetsintrees(est_gts)
    end

    for (i, subset_taxa) in enumerate(subsets)
        output_file = "$(output_file_prefix)_$(i)"
        runtime_file = "$(output_file).runtime"
        output_net_file = "$(output_file).netfile"
        if isfile(output_net_file) && isfile(runtime_file)
            log("SNaQ #$(i)", "Already inferred.")
            continue
        end

        # 1. trim quartets
        subset_q = trim_qt(q, t, subset_taxa)
        df = silently() do
            readTableCF(writeTableCF(subset_q, subset_taxa))
        end

        # 2. infer network (checkpointed)
        nhybrids = retics_in_subnet(true_net, subset_taxa)
        init_tree = pruneTruthFromDecomp(nj_tre, subset_taxa)

        log("SNaQ #$(i)", "Inferring constraint.")
        snaq_runtime = @elapsed snaq_net = silently() do
            # snaq!(init_tree, df, filename=output_file, hmax=nhybrids, seed=seed+i)
            snaq!(init_tree, df, filename=output_file, hmax=nhybrids, seed=seed+i,
                ftolRel=Inf, ftolAbs=Inf, xtolRel=Inf, xtolAbs=Inf, Nfail=3, liktolAbs=1., runs=10)
        end

        # 3. save runtime and network (checkpointed)
        open(runtime_file, "w+") do f
            write(f, "$(snaq_runtime)")
        end
        open(output_net_file, "w+") do f
            write(f, writeTopology(snaq_net))
        end
    end

    # gather up all the constraints and runtimes
    est_constraints::Vector{HybridNetwork} = []
    est_constraint_runtimes::Vector{Float64} = []
    for i=1:length(subsets)
        output_file = "$(output_file_prefix)_$(i)"
        runtime_file = "$(output_file).runtime"
        output_net_file = "$(output_file).netfile"

        push!(est_constraints, readTopology(output_net_file))
        push!(est_constraint_runtimes, parse(Float64, readlines(runtime_file)[1]))
    end

    return est_constraints, est_constraint_runtimes
end


function trim_qt(q, t, subset_taxa)
    temp_q = Vector{typeof(q[1])}([])
    valid_taxonnumbers = Set([i for i=1:length(t) if t[i] in subset_taxa])
    taxonnumber_map = Dict(old_number => findfirst(subset_taxa .== t[old_number]) for old_number in valid_taxonnumbers)
    map_fxn(old_number::Int64) = taxonnumber_map(old_number)

    for quartet in q
        if all(num in valid_taxonnumbers for num in quartet.taxonnumber)
            new_svec = @SVector [taxonnumber_map[quartet.taxonnumber[i]] for i = 1:4]
            quartet_copy = PhyloNetworks.QuartetT{MVector{4, Float64}}(
                length(temp_q) + 1,
                new_svec,
                quartet.data
            )
            push!(temp_q, quartet_copy)
        end
    end

    return temp_q
end
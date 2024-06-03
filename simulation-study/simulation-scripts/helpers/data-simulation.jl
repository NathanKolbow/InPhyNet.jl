# Helper functions for simulating real data
# Uses software in directory `network-merging/software/`
using PhyloCoalSimulations, StaticArraysCore, StaticArrays

function simulate_sequence_data(gts::Vector{HybridNetwork}, output_file_prefix::String, estgt_file::String="", data_dir::String="")
    seq_file_paths = ["$(output_file_prefix)_$(i)" for i=1:length(gts)]
    if all(isfile(f) for f in seq_file_paths)
        @debug "Sequence data files already exist"
        return seq_file_paths
    elseif isfile(estgt_file)
        @debug "Estimated gene trees already exist, skipping sequence generation."
        return seq_file_paths
    end

    # find the correct `-s` flag to use w/ seq-gen
    # s = find_seqgen_s(gts)
    s = 0.036
    

    # simulate sequences w/ seq-gen
    seq_file_paths = run_seqgen_multi(s, gts, output_file_prefix, data_dir = data_dir)
    return seq_file_paths
end


function estimate_gene_trees(seq_files::Vector{String}, estgt_output_file::String)
    if isfile(estgt_output_file)
        @debug "Estimated gene tree file already exists"
        return
    end
    
    estgt_newicks = Array{String}(undef, length(seq_files))
    Threads.@threads for i=1:length(seq_files)
        seq_file = seq_files[i]

        # run iqtree
        run_iqtree(seq_file)

        # read iqtree results
        estgt_newicks[i] = readlines("$(seq_file).treefile")[1]
        clean_est_gt_files("fake_placeholder", seq_file)
    end

    # save iqtree results to file
    open(estgt_output_file, "w+") do f
        for newick in estgt_newicks
            write(f, "$(newick)\n")
        end
    end
end


function silently(f)
    redirect_stdout(Base.DevNull()) do
        redirect_stderr(Base.DevNull()) do 
            return f()
        end
    end
end


"""
Infers constraint networks with taxa in `subsets` given the estimated gene trees in `estgt_file`.
Writes network results to `output_file`, writes runtime results to `output_file.runtimes`
Returns estimated constraint networks and the time to estimate each network.
"""
function infer_constraints(estgt_file::String, output_file::String, subsets::Vector{Vector{String}}, true_net::HybridNetwork; data_dir::AbstractString="/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data/temp_data/")
    if isfile(output_file) && isfile("$(output_file).runtimes")
        constraints = readMultiTopology(output_file)
        runtimes = readlines("$(output_file).runtimes")[1]
        runtimes = [parse(Float64, n) for n in split(runtimes, ",")]
        @debug "Estimated network files already exist"
        return constraints, runtimes
    end
    
    est_gts = readMultiTopology(estgt_file)
    q = t = nothing
    
    redirect_stdout(devnull) do
        q, t = countquartetsintrees(est_gts)
    end

    constraints = Vector{HybridNetwork}([])
    runtimes = Vector{Float64}([])

    for (i, subset_taxa) in enumerate(subsets)
        # reduce (q, t) to only quartets w/ the taxa in `subset`
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

        df = readTableCF(writeTableCF(temp_q, subset_taxa))
        
        # infer w/ snaq
        nhybrids = retics_in_subnet(true_net, subset_taxa)
        temp_out_file = joinpath(data_dir, "snaq_temp_output_$(i)_")
        snaq_net = nothing
        snaq_runtime = -1.
        redirect_stdout(devnull) do
            snaq_runtime = @elapsed snaq_net = snaq!(est_gts[1], df, filename=temp_out_file, hmax=nhybrids, seed=42)
        end
        push!(constraints, snaq_net)
        push!(runtimes, snaq_runtime)
    end

    writeMultiTopology(constraints, output_file)
    open("$(output_file).runtimes", "w+") do f
        for (idx, t) in enumerate(runtimes)
            write(f, idx == length(runtimes) ? "$t" : "$t, ")
        end
    end
    return constraints, runtimes
end


"""
Finds how many reticulations are in the subnetwork of `true_net` composed of the taxa in `subset_taxa`.
Used to determine how many reticulations should be inferred w/ SNaQ.
"""
function retics_in_subnet(true_net, subset_taxa)
    subnet = pruneTruthFromDecomp(true_net, subset_taxa)
    return subnet.numHybrids
end


function load_true_net_ils_adjusted(netid::String, replicatenum::Int64, ils::String)
    truenet = readMultiTopology(getNetworkFilepath(netid))[replicatenum]
    newick = writeTopology(truenet)
    avg_bl = get_avg_bl(truenet)
    newick = "($(newick[1:(length(newick)-1)]):$(avg_bl),OUTGROUP:1.0);"
    truenet = readTopology(newick)

    # ils level : desired average branch length
    # - low       : 2.0
    # - med       : 1.0
    # - high      : 0.5
    # - very high : 0.1

    # 1. make the network have low ILS branch lengths. after adjusting for ILS, we extend branch
    #    lengths such that root to tip distance is unchanged, so this needs to be our baseline
    desired_avg = 2.
    avg_bl = get_avg_bl(truenet)
    for e in truenet.edge
        e.length *= desired_avg / avg_bl
    end

    # 2. record desired final depths
    leaf_depths = get_leaf_depths(truenet)

    # 3. adjust network to desired branch length avg for desired level of ILS
    desired_avg = ils == "low" ? 2. : (ils == "med" ? 1. : (ils == "high" ? 0.5 : 0.1))
    avg_bl = get_avg_bl(truenet)
    for e in truenet.edge
        e.length *= desired_avg / avg_bl
    end
    
    # 4. extend leaves to desired length
    extend_leaves!(truenet, leaf_depths)

    return truenet
end



### Helper functions for the main functions above ###

# Extends the branches above leaves in `net` such that the root to tip distance
# matches that in `desired_depths`.
function extend_leaves!(net::HybridNetwork, desired_depths::AbstractVector{Float64})
    curr_depths = get_leaf_depths(net)
    for (i, adder) in enumerate(desired_depths .- curr_depths)
        edge = getparentedge(net.leaf[i])
        edge.length += adder
    end
end


# Finds the root to tip distance for each tip in net's major tree.
function get_leaf_depths(net::HybridNetwork)
    tre = majorTree(net)
    depths = zeros(net.numTaxa)

    for (i, leaf) in enumerate(net.leaf)
        curr = leaf
        while curr != net.node[net.root]
            depths[i] += getparentedge(curr).length
            curr = getparent(curr)
        end
    end

    return depths
end


function get_descendant_leaves(node::PhyloNetworks.Node)
    q = Vector{PhyloNetworks.Node}([node])
    leaves = Vector{PhyloNetworks.Node}([])
    while length(q) > 0
        curr = pop!(q)
        if curr.leaf
            push!(leaves, curr)
        else
            for child in getchildren(curr)
                push!(q, child)
            end
        end
    end
    return leaves
end


# Finds the `-s` parameter for seq-gen such that gene
# tree estimation error is in a reasonable range
#
# Returns: path to sequence file
function find_seqgen_s(gts::Vector{HybridNetwork}; data_dir::AbstractString="/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data/temp_data/")
    min_s = 0.00001
    max_s = 0.1
    curr_s = 0.003

    desired_gtee = 0.25
    tolerance = 0.01

    ntrees = length(gts)
    gtees = zeros(ntrees)

    while true
        Threads.@threads for i=1:ntrees
            tree = gts[i]
            true_newick = writeTopology(tree)

            temp_seqfile = joinpath(data_dir, "seqgen_$(i).phy")
            temp_gtfile = joinpath(data_dir, "truegt_$(i).treefile")
            writeTopology(tree, temp_gtfile)

            # seq-gen
            run_seqgen(curr_s, temp_gtfile, temp_seqfile)

            # iqtree
            run_iqtree(temp_seqfile)

            # save the result
            est_newick = readlines("$(temp_seqfile).treefile")[1]

            # calculate gtee
            gtee_nrf = calc_gtee(true_newick, est_newick)
            gtees[i] = parse(Float64, gtee_nrf)

            # clean up
            clean_est_gt_files(temp_gtfile, temp_seqfile)
            rm_suppress(temp_seqfile)
        end

        avg_gtee = mean(gtees)
        if avg_gtee - desired_gtee > tolerance
            min_s = curr_s
            curr_s = mean([max_s, min_s])
        elseif avg_gtee - desired_gtee < -tolerance
            max_s = curr_s
            curr_s = mean([max_s, min_s])
        else
            return curr_s
        end
        curr_s = round(curr_s, sigdigits = 4)
    end
end


function run_seqgen(seqgen_s::AbstractFloat, temp_gtfile::AbstractString, temp_seqfile::AbstractString; seq_length::Int64=1000)
    software_path = "/mnt/dv/wid/projects4/SolisLemus-network-merging/software/"
    if Sys.isapple()
        run(pipeline(`$software_path/seq-gen-macos -s$seqgen_s -n1 -f0.3,0.2,0.2,0.3 -mHKY -op -l$(seq_length) $temp_gtfile`, stdout=temp_seqfile, stderr=devnull))
    elseif Sys.islinux()
        run(pipeline(`$software_path/seq-gen-linux -s$seqgen_s -n1 -f0.3,0.2,0.2,0.3 -mHKY -op -l$(seq_length) $temp_gtfile`, stdout=temp_seqfile, stderr=devnull))
    else
        run(pipeline(`$software_path/seq-gen-windows.exe -s$seqgen_s -n1 -f0.3,0.2,0.2,0.3 -t3.0 -mHKY -op -l$(seq_length) $temp_gtfile`, stdout=temp_seqfile, stderr=devnull))
    end
end


# `truegt_file` contains several topologies, and we need to run seq-gen for each of them
function run_seqgen_multi(seqgen_s::AbstractFloat, gts::Vector{HybridNetwork}, output_file_prefix::AbstractString; seq_length::Int64=1000, data_dir::AbstractString="/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data/temp_data/")

    ntrees = length(gts)
    seq_file_paths = Array{String}(undef, ntrees)
    Threads.@threads for i=1:ntrees
        temp_gt_file = joinpath(data_dir, "tempgt_$(i).treefile")
        writeTopology(gts[i], temp_gt_file)

        temp_seq_file = joinpath(data_dir, "tempseq_out_$(i).phy")

        run_seqgen(seqgen_s, temp_gt_file, temp_seq_file, seq_length = seq_length)
        mv(temp_seq_file, "$(output_file_prefix)_$(i)")

        seq_file_paths[i] = "$(output_file_prefix)_$(i)"
        clean_est_gt_files(temp_gt_file, temp_seq_file)
    end
    return seq_file_paths
end


function run_iqtree(temp_seqfile::AbstractString)
    software_path = "/mnt/dv/wid/projects4/SolisLemus-network-merging/software/"
    if Sys.isapple()
        try
            run(pipeline(`$software_path/iqtree-1.6.12-macos -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        catch e
            run(pipeline(`$software_path/iqtree-1.6.12-macos -redo -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        end
    elseif Sys.islinux()
        try
            run(pipeline(`$software_path/iqtree-1.6.12-linux -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        catch e
            run(pipeline(`$software_path/iqtree-1.6.12-linux -redo -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        end
    else
        try
            run(pipeline(`$software_path/iqtree-1.6.12-windows.exe -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        catch e
            run(pipeline(`$software_path/iqtree-1.6.12-windows.exe -redo -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        end
    end
end


function calc_gtee(true_newick::AbstractString, est_newick::AbstractString)
    scriptpath = "/mnt/dv/wid/projects4/SolisLemus-network-merging/software/compare_two_trees.py"
    gtee_nrf = Pipe()
    run(pipeline(`python3 $scriptpath -t1 $true_newick -t2 $est_newick`, stdout=gtee_nrf))
    close(gtee_nrf.in)
    return String(read(gtee_nrf))
end


function clean_iqtree_files(input_seqfile::AbstractString)
    rm_suppress(input_seqfile*".bionj")
    rm_suppress(input_seqfile*".ckp.gz")
    rm_suppress(input_seqfile*".iqtree")
    rm_suppress(input_seqfile*".log")
    rm_suppress(input_seqfile*".mldist")
    rm_suppress(input_seqfile*".model.gz")
    rm_suppress(input_seqfile*".treefile")
end


function clean_est_gt_files(temp_gtfile::AbstractString, temp_seqfile::AbstractString)
    rm_suppress(temp_gtfile)
    rm_suppress(temp_seqfile*".bionj")
    rm_suppress(temp_seqfile*".ckp.gz")
    rm_suppress(temp_seqfile*".iqtree")
    rm_suppress(temp_seqfile*".log")
    rm_suppress(temp_seqfile*".mldist")
    rm_suppress(temp_seqfile*".model.gz")
    rm_suppress(temp_seqfile*".treefile")
end

rm_suppress(file::AbstractString) = try rm(file) catch e end
rm_suppress(files::AbstractArray) = [rmsuppress(file) for file in files]
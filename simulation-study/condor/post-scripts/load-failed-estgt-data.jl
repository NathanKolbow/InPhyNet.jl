include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")
filepath = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/condor/outputs/snaq1_3478_n200r10_2_m15_100_1000_med.out"


function save_all()
    for file in get_valid_file_list()
        save_data(file)
    end
end


function save_data(file_path::String)
    net_id, true_net, rep_num, ngt, ils_level, m, dmethod, seq_len, est_constraints, est_gts, est_constraint_runtimes, est_D, est_namelist = get_data(file_path)

    if estimated_sims_already_performed(net_id, rep_num, ngt, seq_len, ils_level, m)
        printstyled("[ALREADY LOGGED] ", color=:cyan)
        printstyled("$(net_id) #$(rep_num), m=$(m), nloci=$(ngt) seq_len=$(seq_len), ils=$(ils_level)\n", color=:black)
        return
    end

    mnet_time = @elapsed mnet = netnj(est_D, est_constraints, est_namelist)
    if mnet !== nothing && typeof(mnet) <: HybridNetwork
        save_estimated_gts_results(
            net_id, true_net, rep_num, ngt, ils_level, m, dmethod, seq_len,
            mnet, est_constraints, est_gts, est_constraint_runtimes, mnet_time
        )
        printstyled("[SUCCESS] ", color=:green)
        printstyled("$(net_id) #$(rep_num), m=$(m), nloci=$(ngt) seq_len=$(seq_len), ils=$(ils_level)\n", color=:black)
    else
        printstyled("[FAILED] ", color=:red)
        printstyled("$(net_id) #$(rep_num), m=$(m), nloci=$(ngt) seq_len=$(seq_len), ils=$(ils_level)\n", color=:black)
        printstyled("\tnet_id, true_net, rep_num, ngt, ils_level, m, dmethod, seq_len, est_constraints, est_gts, est_constraint_runtimes, est_D, est_namelist=" *
                    "get_data(\"$(file_path)\")")
    end
end


function get_data(file_path::String)
    lines = readlines(file_path)
    idx = 1
    while lines[idx] != "n200r10" && idx <= length(lines)
        idx += 1
    end

    if idx > length(lines)
        return nothing
    end

    net_id = lines[idx]
    true_net = readTopology(lines[idx+1])
    rep_num = parse(Int64, lines[idx+2])
    ngt = parse(Int64, lines[idx+3])
    ils_level = lines[idx+4]
    m = parse(Int64, lines[idx+5])
    dmethod = lines[idx+6]
    seq_len = parse(Int64, lines[idx+7])

    nconstraints = parse(Int64, lines[idx+8])
    est_constraints = Vector{HybridNetwork}([])
    for i=1:nconstraints
        push!(est_constraints, readTopology(lines[idx+8+i]))
    end

    ngts = parse(Int64, lines[idx+8+nconstraints+1])
    est_gts = Vector{HybridNetwork}([])
    for i=1:ngts
        push!(est_gts, readTopology(lines[idx+8+nconstraints+1+i]))
    end

    est_constraint_runtimes = Vector{Float64}([])
    for i=1:nconstraints
        push!(est_constraint_runtimes, parse(Float64, lines[idx+8+nconstraints+1+ngts+i]))
    end
    est_D, est_namelist = calculateAGIC(est_gts)
    est_namelist = Vector{String}(est_namelist)

    return net_id, true_net, rep_num, ngt, ils_level, m, dmethod, seq_len, est_constraints, est_gts, est_constraint_runtimes, est_D, est_namelist
end


function get_valid_file_list()
    base_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/condor/outputs/"
    file_list = Vector{String}([])
    for file in readdir(base_dir, join=true)
        if length(readlines(file)) > 100
            push!(file_list, file)
        end
    end
    return file_list
end





function big_try(netA, netB)
    best_edge = nothing
    best_hwcd = hardwiredClusterDistance(netA, netB, true)

    @info "Starting from $(best_hwcd)"
    for (i, edge) in enumerate(netA.edge)
        if (i % 10) == 0
            print("\t$(i)")
        end

        try
            rootonedge!(netA, edge)
            hwcd = hardwiredClusterDistance(netA, netB, true)
            if hwcd < best_hwcd
                best_edge = edge
                best_hwcd = hwcd
                @info "Improved to $(best_hwcd)"
            end
        catch
        end
    end

    return best_edge, best_hwcd
end

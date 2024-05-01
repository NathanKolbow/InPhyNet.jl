data_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data/"
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")


netfile_names = [filename for filename in readdir(data_dir) if filename[(length(filename)-7):length(filename)] == ".netfile"]
for netfile_name in netfile_names
    netid, replicatenum, ngt, seq_len, ils_level, maxsubsetsize, dmethod = split(netfile_name, "_")[2:8]
    netid = String(netid)
    replicatenum = parse(Int64, replicatenum)
    ngt = parse(Int64, ngt)
    seq_len = parse(Int64, seq_len)
    maxsubsetsize = parse(Int64, maxsubsetsize)
    ils_level = String(ils_level)

    if !estimated_sims_already_performed(netid, replicatenum, ngt, seq_len, ils_level, maxsubsetsize)
        @info "$(netfile_name)"
    end
end

# [ Info: estnets_n200r10_1_1000_500_high_15_AGIC.netfile
est_constraints, est_D, est_namelist = collect_retry_data("estnets_n200r10_1_1000_500_high_15_AGIC.netfile")
mnet = netnj(est_D, est_constraints, est_namelist)
mnet = InPhyNet.netnj_retry_driver(est_D, est_constraints, est_namelist, max_retry_attempts = 1e6)

# [ Info: estnets_n200r10_1_100_500_high_15_AGIC.netfile
est_constraints, est_D, est_namelist = collect_retry_data("estnets_n200r10_1_100_500_high_15_AGIC.netfile")
mnet = netnj(est_D, est_constraints, est_namelist)
mnet = InPhyNet.netnj_retry_driver(est_D, est_constraints, est_namelist, max_retry_attempts = 1e6)

# [ Info: estnets_n200r10_1_100_500_veryhigh_15_AGIC.netfile
est_constraints, est_D, est_namelist = collect_retry_data("estnets_n200r10_1_100_500_veryhigh_15_AGIC.netfile")
mnet = netnj(est_D, est_constraints, est_namelist)
mnet = InPhyNet.netnj_retry_driver(est_D, est_constraints, est_namelist, max_retry_attempts = 1e6)
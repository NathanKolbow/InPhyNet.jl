using PhyloNetworks

# FUNCTIONS FOR LOADING THE DATA EXAMPLES USED IN THE DOCUMENTATION WALKTHROUGH
function load_inphynet_example_gts()::Vector{HybridNetwork}

    path = joinpath(dirname(@__FILE__), "..", "examples", "n50_example_gts.treefile")
    return readmultinewick(path)

end

function loadconflictexample()
	example_dir = joinpath(dirname(pathof(InPhyNet)), "..", "examples", "app5-merge-conflict")
	D_lines = readlines(joinpath(example_dir, "D.csv"))
	namelist = String.(split(D_lines[1], ","))
	D = Matrix{Float64}(undef, 12, 12)
	for i = 1:12
		D[i, :] .= [parse(Float64, value) for value in split(D_lines[i+1], ",")]
	end
	constraints = readmultinewick(joinpath(example_dir, "constraints.net"))
	return D, namelist, constraints
end

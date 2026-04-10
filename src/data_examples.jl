using PhyloNetworks

# FUNCTIONS FOR LOADING THE DATA EXAMPLES USED IN THE DOCUMENTATION WALKTHROUGH
function load_inphynet_example_gts()::Vector{HybridNetwork}

    path = joinpath(dirname(@__FILE__), "..", "examples", "n50_example_gts.treefile")
    return readmultinewick(path)

end


"""
Loads the distance matrix `D`, namelist `namelist`, and constraint networks `constraints`
necessary to reproduce the merge conflict example shown in the Appendices of the InPhyNet
manuscript. Further details are provided upon executing the function in the terminal.
"""
function loadconflictexample()::Nothing
	example_dir = joinpath(dirname(pathof(InPhyNet)), "..", "examples", "app5-merge-conflict")
	D_lines = readlines(joinpath(example_dir, "D.csv"))
	namelist = String.(split(D_lines[1], ","))
	D = Matrix{Float64}(undef, 12, 12)
	for i = 1:12
		D[i, :] .= [parse(Float64, value) for value in split(D_lines[i+1], ",")]
	end
	constraints = readmultinewick(joinpath(example_dir, "constraints.net"))
    Core.eval(Main, :(D = $D))
    Core.eval(Main, :(constraints = $constraints))
    Core.eval(Main, :(namelist = $namelist))

    printstyled("\n\tExample data has been loaded into the variables `D`, `namelist`, and `constraints`. To see the conflict, paste the following into the terminal:\n\n", color=:black)
    printstyled("\t\tinphynet(D, constraints, namelist; refuse_pairwise=true)", color=:cyan)
    printstyled("\n\n\tThis will throw an error stating \"No compatible merge found.\" To use the pairwise version of InPhyNet as a fallback, run the `inphynet` function as usual:\n\n", color=:black)
    printstyled("\t\tinphynet(D, constraints, namelist)", color=:cyan)
    printstyled("\n", color=:black)
    return nothing
end


"""
Loads the distance matrix `D`, namelist `namelist`, and constraint networks `constraints`
necesssary to reproduce the worked example shown in the Appendices of the InPhyNet manuscript.
Further details are provided upon executing the function in the terminal.
"""
function loadworkedexample()::Nothing
	example_dir = joinpath(dirname(pathof(InPhyNet)), "..", "examples", "app4-worked-example")
	D_lines = readlines(joinpath(example_dir, "D.csv"))
	namelist = String.(split(D_lines[1], ","))
	D = Matrix{Float64}(undef, 14, 14)
	for i = 1:14
		D[i, :] .= [parse(Float64, value) for value in split(D_lines[i+1], ",")]
	end
	constraints = readmultinewick(joinpath(example_dir, "constraints.net"))
    Core.eval(Main, :(D = $D))
    Core.eval(Main, :(constraints = $constraints))
    Core.eval(Main, :(namelist = $namelist))

    printstyled("\n\tExample data has been loaded into the variables `D`, `namelist`, and `constraints`. To see the result of the worked example, paste the following into your terminal:\n\n", color=:black)
    printstyled("\t\tinphynet(D, constraints, namelist)", color=:cyan)
    printstyled("\n", color=:black)
    return nothing
end

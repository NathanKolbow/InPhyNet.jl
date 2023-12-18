# Must be run with --project
# args: "<network ID>" "<output filename>" <total sims>

if length(ARGS) != 3
    error("ARGS invalid")
end
netID = ARGS[1]
outfile = abspath(ARGS[2])
totalsims = parse(Int64, ARGS[3])
if !isfile(outfile) error("File $(outfile) not found.") end

include("../pipelines.jl")

truenet, constraints = loadTrueData(netID)
D, namelist = majorinternodedistance(truenet)

# Same as `monophyleticRobustness` from `robustness-fxns.jl` except results are
# written to an output file as they are generated for checkpointing purposes
function monophyleticRobustnessCondor(truenet, constraints, D, namelist, outfilename, uqID; nsim::Int64=25)
    truenewick = writeTopology(truenet)
    
    totalnnigen = Uniform(0, length(constraints) * 4)
    gaussgen = Exponential(1.5)
    out_lk = ReentrantLock()

    df = CSV.read(outfilename, DataFrame)
    if size(df, 1) > 0
        nsim = nsim - sum(df[!,"identifier"] .== uqID)
    end
    if nsim <= 0
        println("All $(sum(df[!,"identifier"] .== uqID)) sims for this uqID already run!")
        return 0
    else
        println("Running $(nsim) sims.")
    end

    Threads.@threads for iter=1:nsim
        # Randomly generate the Gaussian noise parameters
        gaussMean = gaussSd = rand(gaussgen)

        
        # Randomly generate the number of NNI moves
        totalnnimoves = Int64(round(rand(totalnnigen)))
        nnimoves = sample(1:length(constraints), totalnnimoves, replace=true)
        nnimoves = Vector{Int64}([sum(nnimoves .== i) for i=1:length(constraints)])

        esterror = constraintdiffs = nothing
        gausserror = gaussSd
        mnetNewick = ""

        try
            esterror, constraintdiffs, mnetNewick =
                runRobustSim(truenet, constraints, D, namelist, gaussMean, gaussSd, nnimoves)
        catch e
            # findoptQidx having no valid pairs is now caught in `runRobustSim`, so any exceptions
            # here are bad.

            @show typeof(e)
            throw(e)
        end

        lock(out_lk)
        try
            CSV.write(outfilename, DataFrame(
                "identifier" => uqID,
                "merged_newick" => mnetNewick,
                "gauss_sd" => gausserror,
                "total_nni_moves" => totalnnimoves,
                "merged_est_error" => esterror,
                "sum_of_constraint_errors" => sum(constraintdiffs, dims=1)[1,:]
            ), append=true)
        finally
            unlock(out_lk)
        end
    end

    return 0
end

# If we get a 0 return code then we're done!
return monophyleticRobustnessCondor(truenet, constraints, D, namelist, outfile, netID, nsim=totalsims)
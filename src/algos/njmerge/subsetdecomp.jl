# Source for decomposing full species set into subsets
include("mscquartetsinterface.jl")

"""

Helper function that feeds into `decomposeFromQuartets(namelist::ListOfNames, R::QuartetVector)`
"""
function decomposeFromQuartets(filepath::AbstractString; cutoff::Float64=0.01)
    namelist, R = parsequartets(filepath, cutoff=cutoff)
    return decomposeFromQuartets(namelist, R)
end

"""

Decomposes the set of all taxa into disjoint subsets intelligently based on a roughly
estimated species tree and significant quartets from MSCquartets.

## Arguments

- T: the quickly estimated, rough species tree
- R: the set of statistically significant quartets, as obtained by MSCquartets

## Return Value

Returns a tuple. The first entry is list of subsets that are required in order to indentify
all relevant reticulations. The second entry is the remaining taxa that aren't in the first entry.
"""
function decomposeFromQuartets(namelist::ListOfNames, R::QuartetVector)
    check_quartets(R)
    reqsubsets = requiredhybridsubsets(namelist, R)
    treetaxa = [name for name in namelist if !inany(name, reqsubsets)]
    return (reqsubsets, treetaxa)
end


"""

Helper function that finds whether `name` is listed in any subset in `subsets`
"""
function inany(name::AbstractString, subsets::QuartetVector)
    for subset in subsets
        if name in subset return true end
    end
    return false
end


"""

Current base subset decomposition schema. First, we find the smallest subset `S` of names such that
all quartets have at least one taxa in the subset. Then, for each entry of `S`, we gather the set of
3 other taxa it appears with most often in the quartets in `R`.

The idea behind these subsets is to guarantee identifiability of hybridizations, hence why they are
required.
"""
function requiredhybridsubsets(namelist::ListOfNames, R::QuartetVector; cutoff::Float64=0.01)
    minsubset = minimumuniquesubset(namelist, R, cutoff=cutoff)
    reqsubsets = kconnected(namelist, minsubset, R, k=4, cutoff=0.01)

    i = 1
    j = 2
    while i < length(reqsubsets)
        un = union(reqsubsets[i], reqsubsets[j])
        
        if length(un) != length(reqsubsets[i]) + length(reqsubsets[j])
            reqsubsets[i] = un
            deleteat!(reqsubsets, j)
            j = i + 1
        else
            j += 1
        end

        if j > length(reqsubsets)
            i += 1
            j = i + 1
        end
    end

    return reqsubsets
end


"""

Finds the smallest subset `S` of names in `namelist` such that every quartet in `R` has at least
one of its taxa in `S`.
"""
function minimumuniquesubset(namelist::ListOfNames, R::QuartetVector; cutoff::Float64=0.01)
    subset = Vector{AbstractString}()

    while length(R) > 0
        inquartet = [[name in quartet for quartet in R] for name in namelist]
        counts = sum.(inquartet)
        addidx = argmax(counts)
        addname = namelist[addidx]

        rowfilter = filter(i -> i != addidx, 1:length(namelist))
        namelist = view(namelist, rowfilter)
        R = view(R, .!inquartet[addidx])
        push!(subset, addname)
    end

    return subset
end


"""

Takes `minsubset`, a subset defined in the function `minimumuniquesubset`, and, for each entry `i` of `minsubset`,
gathers the `k-1` taxa that `i` appears with most often across all quartets.
"""
function kconnected(namelist::ListOfNames, minsubset::Vector{<:AbstractString}, R::QuartetVector; k::Int64=4, cutoff=0.01)
    k > 1 || error("k must be > 1")
    kconn = Vector{Vector{String}}()

    for name in minsubset
        tempR = view(R, [name in quartet for quartet in R])
        counts = [sum([partnername in quartet for quartet in tempR]) for partnername in namelist]
        topk = sortperm(counts, rev=true)[1:min(k, length(counts))]
        push!(kconn, namelist[topk])
    end

    return kconn
end


"""

Makes sure the vector of abstract strings (that we call quartets) each have 4 unique entries.
"""
function check_quartets(R::QuartetVector)
    for quartet in R
        if length(quartet) != 4
            throw(ArgumentError("Each quartet must have exactly 4 entries."))
        elseif length(unique(quartet)) != 4
            throw(ArgumentError("Quartets must have unique entries."))
        end
    end
end


################################################################
# These functions are actually related to subset decomposition #
################################################################

function requiredsubsets(df, cutoff)
    minsubset = minimumuniquesubset(df, cutoff=cutoff)
    return kconnected(df, minsubset, k=4, cutoff=0.01)
end

function minimumuniquesubset(df; cutoff=0.01)
    df = deepcopy(df)
    N = sum(df[:,"HBp_T3"] .<= cutoff)
    subset = []
    namelist = names(df)[2:(findfirst(names(df) .== "12|34")-1)]

    while sum(df[:,"HBp_T3"] .<= cutoff) > 0
        reject = df[:,"HBp_T3"] .<= cutoff
        counts = [sum(df[:,name] .== 1 .&& reject) for name in namelist if name in names(df)]
        
        removeidx = argmax(counts)
        myname = namelist[removeidx]

        numremoved = sum(df[:,myname] .== 1 .&& reject)
        removerows = df[:,myname] .== 1
        df = df[.!removerows,:]

        push!(subset, namelist[removeidx])
        deleteat!(namelist, removeidx)
        deleteat!(counts, removeidx)
    end

    return subset
end

function kconnected(df, minsubset; k::Int64=4, cutoff=0.01)
    k > 1 && k <= 4 || error("k must be in [2, 4]")
    kconnected = []
    namelist = names(df)
    namelist = namelist[2:(findfirst(namelist .== "12|34")-1)]

    for name in minsubset
        tempdf = view(df, df[:,name] .== 1 .&& df[:,"HBp_T3"] .<= cutoff, :)
        counts = [sum(tempdf[:,partnername] .== 1) for partnername in namelist]
        topk = sortperm(counts, rev=true)[1:k]
        push!(kconnected, namelist[topk])
    end
    return kconnected
end

################################################################
#                                                              #
################################################################
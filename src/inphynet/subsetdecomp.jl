# Source for decomposing full species set into subsets

"""

Helper function that feeds into `decomposeFromQuartets(namelist::ListOfNames, R::QuartetVector)`
"""
function decomposeFromQuartets(filepath::AbstractString; cutoff::Float64=0.01, kwargs...)
    if !isfile(filepath)
        throw(ArgumentError("File $filepath not found."))
    elseif cutoff < 0 || cutoff > 1
        throw(ArgumentError("Cutoff must be on range [0, 1]."))
    end

    namelist, R = parsequartets(filepath, cutoff=cutoff)
    return decomposeFromQuartets(namelist, R; kwargs...)
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
function decomposeFromQuartets(namelist::ListOfNames, R::QuartetVector; maxexpansionlen::Int64=9,
    distmat::Union{Matrix{Float64},Nothing}=nothing)
    check_quartets(R)
    hybsubsets = requiredhybridsubsets(namelist, R)
    treetaxa = [name for name in namelist if !inany(name, hybsubsets)]

    if distmat !== nothing
        expandhybsubsetsD!(hybsubsets, treetaxa, R, distmat, namelist, maxlen=maxexpansionlen)
    else
        expandhybsubsets!(hybsubsets, treetaxa, R, maxlen=maxexpansionlen)
    end

    # hyblengths = [length(sub) for sub in hybsubsets]
    # while minimum(hyblengths) < minexpansionlen
    #     minidx1 = argmin(hyblengths)
    #     minidx2 = argmin(hyblengths[filter(i -> i != minidx1, 1:length(hybsubsets))])
    #     hybsubsets[minidx1] = [hybsubsets[minidx1]; hybsubsets[minidx2]]
    #     deleteat!(hybsubsets, minidx2)
    #     hyblengths = [length(sub) for sub in hybsubsets]
    # end

    return (hybsubsets, treetaxa)
end


"""
    expandhybsubsetsD!(hybsubs, other, quartets, D; maxlen::Int64=9)

Helper function; `hybsubs` is the subset of taxa selected for hybridization identification,
`other` is all other taxa which are currently unplaced. This function expands the taxa in
`other` appropriately into various subsets in `hybsubs` based on distance matrix `D`.
"""
function expandhybsubsetsD!(hybsubs, others, quartets, D, Dnames; maxlen::Int64=9, minlen::Int64=6)
    removeMaxLenSets(hsubs) = view(hsubs, filter(i -> length(hsubs[i]) < maxlen, 1:length(hsubs)))
    hybsubs = removeMaxLenSets(hybsubs)
    nametoidx(name) = findfirst(name .== Dnames)

    # If at the end of the `while` loop `others` is non-zero and less than `minlen`,
    # we need to undo our last few placements so that everything has length
    # greater than minlen, so we log our placements
    placements = []

    while length(hybsubs) != 0 && length(others) != 0
        otheridx = nametoidx(others[1])
        hybidxs = [[nametoidx(name) for name in subset] for subset in hybsubs]
        avgdists = [mean(D[idxs, otheridx]) for idxs in hybidxs]
        moveto = argmin(avgdists)
        push!(hybsubs[moveto], others[1])
        deleteat!(others, 1)
        push!(placements, (hybsubs[moveto], length(hybsubs[moveto])))

        hybsubs = removeMaxLenSets(hybsubs)
    end

    if length(others) != 0 && length(others) < minlen
        i = length(placements)
        while length(others) < minlen
            if i == 0
                @warn "Could not find suitable subset decomposition for the `maxexpansionlen` and `minexpansionlen` arguments provided."
                return
            end

            record = placements[i]
            if length(record[1]) > minlen
                push!(others, record[1][record[2]])
                deleteat!(record[1], record[2])
            end

            i -= 1
        end
    end
end


function expandhybsubsets!(hybsubs, others, quartets; maxlen::Int64=9)
    hybsubs = view(hybsubs, filter(i -> length(hybsubs[i]) < maxlen, 1:length(hybsubs)))
    
    movedfromothers = []
    for (i, othertaxa) in enumerate(others)
        if length(hybsubs) == 0 break end

        counts = []
        for hybsubset in hybsubs
            total = 0
            for q in quartets
                if any([hybtax in q for hybtax in hybsubset])
                    total += 1
                end
            end
            push!(counts, total)
        end
        maxcount = argmax(counts)
        push!(hybsubs[maxcount], othertaxa)
        push!(movedfromothers, othertaxa)

        if length(hybsubs[maxcount]) >= maxlen
            hybsubs = view(hybsubs, filter(i -> i != maxcount, 1:length(hybsubs)))
        end
    end
    for del in movedfromothers
        deleteat!(others, findfirst(others .== del))
    end
end


"""

Helper function that finds whether `name` is listed in any subset in `subsets`
"""
@inline function inany(name::AbstractString, subsets::QuartetVector)
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
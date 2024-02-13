"""

Reads output from MSCquartets in R and parses it so that it can be used efficiently with our julia code.
"""
function parsequartets(filepath::AbstractString; cutoff::Float64=0.01)
    return parsequartets(CSV.read(filepath, DataFrame), cutoff=cutoff)
end

"""

Reads output from MSCquartets in R and parses it so that it can be used efficiently with our julia code.
"""
function parsequartets(ptable::DataFrame; cutoff::Float64=0.01)
    check_ptable(ptable)

    if names(ptable)[1] == "Column1"
        ptable = view(ptable, :, 2:size(ptable, 2))
    end

    namelist = names(ptable)
    pvalcol = findfirst(namelist .== "HBp_T3")
    namecols = 1:(findfirst(namelist .== "12|34")-1)
    namelist = namelist[namecols]
    ptable = Matrix(ptable)
    qlist = QuartetVector()
    
    for row in 1:size(ptable, 1)
        if ptable[row, pvalcol] <= cutoff
            taxaidxs = findall(ptable[row, namecols] .== 1)
            push!(qlist, Quartet(namelist[taxaidxs]))
        end
    end

    return namelist, qlist
end

"""

Checks the input table to make sure that it is Holmes-Bonferroni adjusted.
TODO: if it isn't, manually adjust it instead of throwing an error.
"""
check_ptable(ptable::DataFrame) = 
    ("HBp_T3" in names(ptable)) || throw(ArgumentError("MSCquartets p-table must be Holmes-Bonferroni adjusted."))

# parsequartets("C:\\Users\\Nathan\\repos\\network-merging\\src\\tests\\netnjmerge\\10000-constraints1.hbptab")
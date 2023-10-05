# Calculates the probability of coalescent outcome `sets`
# when there are `N` total lineages in a branch of length `bl`
# function _calculatecoalescentprobability(sets, N::Real, bl::Real)
#     if length(sets) == N
#         return exp(-N*bl)
#     end
#     return 0
# end

# TODO: instead of setting as globals, keep these variables in the main algo and pass them to this function
# TODO: make sure all inputs to `exp` functions and all else are `BigFloat`s
global const savedprobs = Ref{Dict{Tuple{Real, Real, Real}, BigFloat}}(Dict{Tuple{Real, Real, Real}, BigFloat}())
global const savedcomps = Ref{Int64}(0)
"""
    _calculatecoalescentprobability(N::Real, O::Real, bl::Real)

Calculate the probability of `N` lineages coalescing into `O` lineages
in a branch of length `bl`. Real values are accepted for `N` and `O` but
they must be castable to integer.

# Arguments
- `N::Real`: the number of "i**N**put" lineages; the nubmer of lineages at the beginning of the branch
- `O::Real`: the number of "**O**utput" lineages; the nubmer of lineages at the beginning of the branch
- `bl::Real`: length of the branch that the coalescence is happening in

# Return Value
Real value in (0, 1]; the probability of the coalescence under the given parameters happening.
"""
function _calculatecoalescentprobability(N::Real, O::Real, bl::Real)
    # println("("*string(N)*","*string(O)*","*string(bl)*")")
    # Make sure `N` and `O` are castable as `Real`
    try 
        N = Int64(N)
        O = Int64(O)
        bl = BigFloat(bl)
    catch err
        if typeof(err) <: MethodError
            error("Either `N` or `O` are not convertible to `Int64` in method call: _calculatecoalescentprobability("*string(N)*", "*string(O)*", "*string(bl)*")")
        else
            rethrow(err)
        end
    end

    N >= O > 0 || error("`N` must be >= `O` and both must be greater than 0")
    savedict = savedprobs[]

    # if haskey(savedict, (N, O, bl))
    #     savedcomps[] += 1
    #     return savedict[N, O, bl]
    # end

    # Pre-computed probas
    # functions for all outcomes with N=4: https://www.desmos.com/calculator/0x4jmv0muv
    if N == O
        savedict[N, O, bl] = exp(-binomial(N,2)*bl)
        return exp(-binomial(N,2)*bl)
    elseif N - O == 1
        c1 = binomial(N, 2)
        c2 = binomial(O, 2)
        savedict[N, O, bl] = (exp(-bl*c1) - exp(-bl*c2)) / (c2 - c1)

        return (exp(-bl*c1) - exp(-bl*c2)) / (c2 - c1)
    elseif N - O == 2
        c1 = binomial(N, 2)
        c2 = binomial(N-1, 2)
        c3 = binomial(O, 2)
        
        # TODO: simplify the computation below
        num = ((c2-c1)*exp(c2*bl+c1*bl)+((c3-c2)*exp(c2*bl)+(c1-c3)*exp(c1*bl))*exp(c3*bl))*exp(-(c1+c2+c3)*bl)
        denom = (c2-c1)*(c3^2-(c1+c2)*c3+c1*c2)
        savedict[N, O, bl] = num/denom

        return num/denom
    elseif N - O == 3
        c1 = binomial(N, 2)
        c2 = binomial(N-1, 2)
        c3 = binomial(N-2, 2)
        c4 = binomial(O, 2)

        # TODO: this computation is wrong
        num = (((((exp(bl*c2)-exp(bl*c1))*c3-c2*exp(bl*c2)+c1*exp(bl*c1))*exp(bl*c3)+(exp(bl*c1)*c2-c1*exp(bl*c1))*exp(bl*c2))*c4^2+(((exp(bl*c1)-exp(bl*c2))*c3^2+c2^2*exp(bl*c2)-c1^2*exp(bl*c1))*exp(bl*c3)+(c1^2*exp(bl*c1)-exp(bl*c1)*c2^2)*exp(bl*c2))*c4+((c2*exp(bl*c2)-c1*exp(bl*c1))*c3^2+(c1^2*exp(bl*c1)-c2^2*exp(bl*c2))*c3)*exp(bl*c3)+(c1*exp(bl*c1)*c2^2-c1^2*exp(bl*c1)*c2)*exp(bl*c2))*exp(bl*c4)+((c1*exp(bl*c1)-exp(bl*c1)*c2)*exp(bl*c2)*c3^2+(exp(bl*c1)*c2^2-c1^2*exp(bl*c1))*exp(bl*c2)*c3+(c1^2*exp(bl*c1)*c2-c1*exp(bl*c1)*c2^2)*exp(bl*c2))*exp(bl*c3))*exp(-(c1+c2+c3+c4)*bl)
        denom = (c2-c1)*(c3^2-(c1+c2)*c3+c1*c2)*(c4^3-(c1+c2+c3)*c4^2+((c1+c2)*c3+c1*c2)*c4-c1*c2*c3)

        savedict[N, O, bl] = num/denom
        return num/denom
    end

    error("Coalescent probabilities only implemented up to an in-out difference of 3; received (N,O) := ("*string(N)*","*string(O)*")")
end
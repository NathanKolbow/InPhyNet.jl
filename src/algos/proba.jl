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
function _calculatecoalescentprobability(N::Real, O::Real, bl::Real; complog::Union{CompDict,Nothing}=nothing)
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

    if complog !== nothing && haskey(complog, (N, O, bl))
        return complog[N, O, bl]
    end

    # Pre-computed probas
    # functions for all outcomes with N=4: https://www.desmos.com/calculator/0x4jmv0muv
    retval = 0.
    if N == O
        retval = exp(-binomial(N,2)*bl)
    elseif N - O == 1
        c1 = binomial(N, 2)
        c2 = binomial(O, 2)

        retval = (exp(-bl*c1) - exp(-bl*c2)) / (c2 - c1)
    elseif N - O == 2
        c1 = binomial(N, 2)
        c2 = binomial(N-1, 2)
        c3 = binomial(O, 2)
        
        c1bl = exp(c1*bl)
        c2bl = exp(c2*bl)
        c3bl = exp(c3*bl)
        
        # num = ((c2-c1)*exp(c2*bl+c1*bl)+((c3-c2)*exp(c2*bl)+(c1-c3)*exp(c1*bl))*exp(c3*bl))*exp(-(c1+c2+c3)*bl)
        num = ((c2-c1)*c1bl*c2bl+((c3-c2)*c2bl+(c1-c3)*c1bl)*c3bl)
        denom = (c2-c1)*(c3^2-(c1+c2)*c3+c1*c2)*c1bl*c2bl*c3bl

        retval = num/denom
    elseif N - O == 3
        c1 = binomial(N, 2)
        c2 = binomial(N-1, 2)
        c3 = binomial(N-2, 2)
        c4 = binomial(O, 2)

        c1bl = exp(c1*bl)
        c2bl = exp(c2*bl)
        c3bl = exp(c3*bl)
        c4bl = exp(c4*bl)

        # num = (((((exp(bl*c2)-exp(bl*c1))*c3-c2*exp(bl*c2)+c1*exp(bl*c1))*exp(bl*c3)+(exp(bl*c1)*c2-c1*exp(bl*c1))*exp(bl*c2))*c4^2+(((exp(bl*c1)-exp(bl*c2))*c3^2+c2^2*exp(bl*c2)-c1^2*exp(bl*c1))*exp(bl*c3)+(c1^2*exp(bl*c1)-exp(bl*c1)*c2^2)*exp(bl*c2))*c4+((c2*exp(bl*c2)-c1*exp(bl*c1))*c3^2+(c1^2*exp(bl*c1)-c2^2*exp(bl*c2))*c3)*exp(bl*c3)+(c1*exp(bl*c1)*c2^2-c1^2*exp(bl*c1)*c2)*exp(bl*c2))*exp(bl*c4)+((c1*exp(bl*c1)-exp(bl*c1)*c2)*exp(bl*c2)*c3^2+(exp(bl*c1)*c2^2-c1^2*exp(bl*c1))*exp(bl*c2)*c3+(c1^2*exp(bl*c1)*c2-c1*exp(bl*c1)*c2^2)*exp(bl*c2))*exp(bl*c3))*exp(-(c1+c2+c3+c4)*bl)
        num = (((((c2bl-c1bl)*c3-c2*c2bl+c1*c1bl)*c3bl+(c1bl*c2-c1*c1bl)*c2bl)*c4^2+(((c1bl-c2bl)*c3^2+c2^2*c2bl-c1^2*c1bl)*c3bl+(c1^2*c1bl-c1bl*c2^2)*c2bl)*c4+((c2*c2bl-c1*c1bl)*c3^2+(c1^2*c1bl-c2^2*c2bl)*c3)*c3bl+(c1*c1bl*c2^2-c1^2*c1bl*c2)*c2bl)*c4bl+((c1*c1bl-c1bl*c2)*c2bl*c3^2+(c1bl*c2^2-c1^2*c1bl)*c2bl*c3+(c1^2*c1bl*c2-c1*c1bl*c2^2)*c2bl)*c3bl)
        denom = (c2-c1)*(c3^2-(c1+c2)*c3+c1*c2)*(c4^3-(c1+c2+c3)*c4^2+((c1+c2)*c3+c1*c2)*c4-c1*c2*c3)*c1bl*c2bl*c3bl*c4bl

        retval = num/denom
    else
        # Pass the job to the less efficient algorithm 
        # that can work on arbitrary (N,O)
        retval = g(N, O, bl)
    end

    if complog !== nothing complog[N, O, bl] = retval end
    return retval
end


# TODO: push `cs` into this so we don't redo a ton of binomials
function _calculatetotalcoalescentprobability(N::Real, O::Real, bl::Real; complog::Union{CompDict,Nothing}=nothing, cs::Nothing=nothing)
    return _calculatecoalescentprobability(N, O, bl, complog=complog) * prod([binomial(i, 2) for i=(O+1):N])
end
# Calculates the probability of coalescent outcome `sets`
# when there are `N` total lineages in a branch of length `bl`
# function _calculatecoalescentprobability(sets, N::Real, bl::Real)
#     if length(sets) == N
#         return exp(-N*bl)
#     end
#     return 0
# end

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
    # Make sure `N` and `O` are castable as `Real`
    try 
        N = Int64(N)
        O = Int64(O)
    catch err
        if typeof(err) <: MethodError
            error("Either `N` or `O` are not convertible to `Int64` in method call: _calculatecoalescentprobability("*string(N)*", "*string(O)*", "*string(bl)*")")
        else
            rethrow(err)
        end
    end

    N >= O > 0 || error("`N` must be >= `O` and both must be greater than 0")

    # Pre-computed probas
    # functions for all outcomes with N=3: https://www.desmos.com/calculator/atrnmeryug
    if N == O
        return exp(-binomial(N,2)*bl)
    elseif N - O == 1
        c1 = binomial(N, 2)
        c2 = binomial(O, 2)

        return (exp(-bl*c1) - exp(-bl*c2)) / (c2 - c1)
    elseif N - O == 2
        c1 = binomial(N, 2)
        c2 = binomial(N-1, 2)
        c3 = binomial(O, 2)

        # TODO: simplify the computation below
        return ((c2-c1)*exp(c2*bl+c1*bl)+((c3-c2)*exp(c2*bl)+(c1-c3)*exp(c1*bl))*exp(c3*bl))*exp(-(c1+c2+c3)*bl) / ((c2-c1)*(c3^2-(c1+c2)*c3+c1*c2))
        # return (c1*exp(c2*bl+c1*bl)+(((c3-c2)*exp(c1*bl)-c3+c2)*exp(c2*bl)-c1*exp(c1*bl))*exp(c3*bl))*exp(-(c1+c2+c3)*bl) / (c1*(c2-c1)*(c3-c2))

        # Only correct for N=3
        # return (((exp(bl*c1)-1)*c2-c1*exp(bl*c1))*exp(bl*c2)+c1*exp(bl*c1))*exp(-bl*c2-bl*c1) / (c1*c2*(c2-c1))
    elseif N - O == 3
        c1 = binomial(N, 2)
        c2 = binomial(N-1, 2)
        c3 = binomial(N-2, 2)
        c4 = binomial(O, 2)

        num = (((((exp(bl*c2)-exp(bl*c1))*c3-c2*exp(bl*c2)+c1*exp(bl*c1))*exp(bl*c3)+(exp(bl*c1)*c2-c1*exp(bl*c1))*exp(bl*c2))*c4^2+(((exp(bl*c1)-exp(bl*c2))*c3^2+c2^2*exp(bl*c2)-c1^2*exp(bl*c1))*exp(bl*c3)+(c1^2*exp(bl*c1)-exp(bl*c1)*c2^2)*exp(bl*c2))*c4+((c2*exp(bl*c2)-c1*exp(bl*c1))*c3^2+(c1^2*exp(bl*c1)-c2^2*exp(bl*c2))*c3)*exp(bl*c3)+(c1*exp(bl*c1)*c2^2-c1^2*exp(bl*c1)*c2)*exp(bl*c2))*exp(bl*c4)+((c1*exp(bl*c1)-exp(bl*c1)*c2)*exp(bl*c2)*c3^2+(exp(bl*c1)*c2^2-c1^2*exp(bl*c1))*exp(bl*c2)*c3+(c1^2*exp(bl*c1)*c2-c1*exp(bl*c1)*c2^2)*exp(bl*c2))*exp(bl*c3))*exp(-(c1+c2+c3+c4)*bl)
        denom = (c2-c1)*(c3^2-(c1+c2)*c3+c1*c2)*(c4^3-(c1+c2+c3)*c4^2+((c1+c2)*c3+c1*c2)*c4-c1*c2*c3)

        return num/denom
    end

    error("Coalescent probabilities only implemented up to an in-out difference of 3; received (N,O) := ("*string(N)*","*string(O)*")")
end

_calculatecoalescentprobability(3, 3, 1.) + 3*_calculatecoalescentprobability(3, 2, 1.) + 3*_calculatecoalescentprobability(3, 1, 1.)

_calculatecoalescentprobability(4, 4, 1.) + 
    binomial(4, 2)*_calculatecoalescentprobability(4, 3, 1.) +
    binomial(4, 2)*binomial(4-1,2)*_calculatecoalescentprobability(4, 2, 1.) +
    binomial(4, 2)*binomial(4-1,2)*binomial(4-2,2)*_calculatecoalescentprobability(4, 1, 1.)
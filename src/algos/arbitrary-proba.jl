# Code for calculating coal probas w/ arbitrary (N, O)

# TODO: probably just replace `a` and `b` with `exp` for the exponent
#       if we're going to be storing the evaluated `c` values then we
#       may as well just use 1 field rather than 2
# TODO: don't track `t` in the structs, that keeps track of itself
#       in the main for loop
struct foneobj      # e.g. 1/(c4-c2) [f(t5, c2, c4) - 1]
    t::Int64        # e.g.               5
    a::BigInt       # e.g.                  c2 (full value of c2)
    b::BigInt       # e.g.                      c4 (full value of c4)
    invcoef::BigInt # e.g. 1/(c4-c2), stored as (c4-c2)
end
sign(f::foneobj) = ifelse(f.invcoef > 0, 1, -1)
isneg(f::foneobj) = f.invcoef < 0

Base.:-(f::foneobj) = foneobj(f.t, f.a, f.b, -f.invcoef)
function Base.:*(f1::foneobj, f2::foneobj)
    f1.t == f2.t || error("f1 and f2 have different t's")
    
    # Some if/elses to help save big multiplications
    # TODO: benchmark whether this actually saves any time
    firstprod = nothing
    if f1.a == f2.b
        firstprod = fobj(f1.t, f1.b, f2.a, f1.invcoef * f2.invcoef)
    elseif f2.b == f1.a
        firstprod = fobj(f1.t, f2.a, f1.b, f1.invcoef * f2.invcoef)
    else
        firstprod = fobj(f1.t, f1.a + f2.a, f1.b + f2.b, f1.invcoef * f2.invcoef)
    end

    return [
        firstprod,
        fobj(f1.t, f1.a, f1.b, -f1.invcoef),
        fobj(f2.t, f2.a, f2.b, -f2.invcoef),
        1.
    ]
end
Base.:*(f::foneobj, r::Real) = foneobj(f.t, f.a, f.b, r*f.invcoef)
Base.:*(r::Real, f::foneobj) = f*r


struct fobj
    t::Int64
    a::BigInt
    b::BigInt
    invcoef::BigInt
end
function Base.:*(f1::foneobj, f::fobj)
    f1.t == f.t || error("f1 and f have different t's")
    
    # Some if/elses to help save big multiplications
    # TODO: benchmark whether this actually saves any time
    firstprod = nothing
    if f1.a == f.b
        firstprod = fobj(f1.t, f1.b, f.a, f1.invcoef * f.invcoef)
    elseif f.b == f1.a
        firstprod = fobj(f1.t, f.a, f1.b, f1.invcoef * f.invcoef)
    else
        firstprod = fobj(f1.t, f1.a + f.a, f1.b + f.b, f1.invcoef * f.invcoef)
    end

    return [
        firstprod,
        fobj(f.t, f.a, f.b, -f1.invcoef)
    ]
end
Base.:*(f1::fobj, f2::fobj) = fobj(f1.t, f1.a+f2.a, f1.b+f2.b, f1.invcoef*f2.invcoef)
Base.:*(f::fobj, f1::foneobj) = f1*f
Base.:-(f::fobj) = fobj(f.t, f.a, f.b, -f.invcoef)



#############################
# All integration functions #
#############################

function g(N::Integer, O::Integer, bl::BigFloat)
    N - O > 1 || error("N - O > 1 required.")

    # Pre-compute necessary c's
    cs = Array{BigInt}(undef, N-O+1)
    for i=1:(N-O+1) cs[i] = binomial(N-i+1, 2) end
    
    prevg = foneobj(2, cs[1], cs[2], cs[2]-cs[1])
    for i=2:(N-O)
        # At each iter we have to solve `int_0^{t(i+1)} g(i-1)f(ti, ci, c(i+1)) dti`
        # g(i-1) is a list of foneobjs, stored in prevres
        prevg = coalintegrate(prevg, i, cs)
    end
    
    return sum([evaluate(f, bl) for f in prevg]) * exp(-cs[N-O+1]*bl)
end
g(N::Integer, O::Integer, bl::Real) = g(N, O, BigFloat(bl))

coalintegrate(g::Union{fobj, foneobj}, iter::Integer, cs::Array{BigInt}) = coalintegrate([g], iter, cs)
function coalintegrate(g::AbstractArray{<:Union{fobj, foneobj}}, iter::Integer, cs::Array{BigInt})
    # Solving: int_0^{t_{iter+1}} g(iter-1)f(t_{iter}, c_{iter}, c_{iter+1})
    
    # Expand into multiple integrals
    fones = Vector{foneobj}()

    iterf = fobj(iter, cs[iter], cs[iter+1], 1)
    for (i, fone) in enumerate(g)
        products = fone * iterf
        products = ifelse(typeof(products) <: AbstractArray, products, [products])

        for (j, fprod) in enumerate(products)
            push!(fones, fintegrate(fprod))
        end
    end

    return fones
end
fintegrate(f::fobj) = foneobj(f.t+1, f.a, f.b, f.invcoef * (f.b - f.a))
evaluate(f::foneobj, bl::BigFloat) = BigFloat(1.) / f.invcoef * (exp((f.b - f.a) * bl) - 1)
evaluate(f::foneobj, bl::Real) = evaluate(f, BigFloat(bl))





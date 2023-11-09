import Base: getindex, setindex!

const EdgeOrNA = Union{Edge, Nothing}
struct ReticMap
    map::Dict{Edge, Tuple{EdgeOrNA, EdgeOrNA}}
    function ReticMap(constraints::Vector{HybridNetwork})
        # map all existing reticulations in `constraints` to `(nothing, nothing)`
        d = Dict()
        for net in constraints
            for e in net.edge
                if e.hybrid && !e.isMajor
                    d[e] = (nothing, nothing)
                end
            end
        end
        new(d)
    end
end

Base.getindex(r::ReticMap, e::Edge) = r[e]
Base.getindex(r::ReticMap, e::Edge, dir::AbstractString) = ifelse(dir == "from", r[e][1], r[e][2])
function Base.setindex!(r::ReticMap, e::Edge, dir::AbstractString, newe::Edge)
    if dir == "fron"
        return r[e] = (newe, r[e][2])
    else
        return r[e] = (r[e][1], newe)
    end
end
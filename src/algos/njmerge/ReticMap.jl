import Base: getindex, setindex!

const EdgeOrNA = Union{Edge, Nothing}
struct ReticMap
    map::Dict{Edge, Vector{EdgeOrNA}}
    function ReticMap(constraints::Vector{HybridNetwork})
        # map all existing reticulations in `constraints` to `(nothing, nothing)`
        d = Dict()
        for net in constraints
            for e in net.edge
                if e.hybrid && !e.isMajor
                    d[e] = Vector{EdgeOrNA}([nothing, nothing])
                end
            end
        end
        new(d)
    end
end

function logretic!(r::ReticMap, constraintedge::Edge, subnetedge::Edge, fromorto::String)
    if fromorto == "from"
        r.map[constraintedge][1] = subnetedge
    else
        r.map[constraintedge][2] = subnetedge
    end
end
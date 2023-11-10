import Base: getindex, setindex!

struct ReticMap
    map::Dict{Edge, Vector{Edge}}
    function ReticMap(constraints::Vector{HybridNetwork})
        # map all existing reticulations in `constraints` to `(nothing, nothing)`
        d = Dict()
        for net in constraints
            for e in net.edge
                if e.hybrid && !e.isMajor
                    d[e] = Vector{Edge}()
                end
            end
        end
        new(d)
    end
end

logretic!(r::ReticMap, constraintedge::Edge, subnetedge::Edge) = push!(r.map[constraintedge], subnetedge)
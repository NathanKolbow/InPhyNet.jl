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
        if !(r.map[constraintedge][1] === nothing) throw(ErrorException("Overriding from edge")) end
        r.map[constraintedge][1] = subnetedge
        if r.map[constraintedge][2] == subnetedge
            @error("equiv edges A")
            throw(ErrorException("equiv edges A"))
        end
    else
        if !(r.map[constraintedge][2] === nothing) throw(ErrorException("Overriding to edge")) end
        r.map[constraintedge][2] = subnetedge
        if r.map[constraintedge][1] == subnetedge
            @error("equiv edges B")
            throw(ErrorException("equiv edges B"))
        end
    end
end

function check_reticmap(r::ReticMap)
    for (i, key) in enumerate(keys(r.map))
        if length(r.map[key]) != 2
            error("ReticMap key $i has $(length(r.map[key])) attached edges.")
        elseif sum(r.map[key] .!== nothing) != 2
            println(r.map[key])
            error("ReticMap key $i has $(sum(r.map[key] .!== nothing)) attached non-nothing edges.")
        end
    end
end
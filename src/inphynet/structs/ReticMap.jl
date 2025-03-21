const EdgeOrNA = Union{Edge, Nothing}
struct ReticMap
    map::Dict{Node, Vector{EdgeOrNA}}
    edge_map::Dict{Edge, Node}

    function ReticMap(constraints::Vector{HybridNetwork})
        # map all existing reticulations in `constraints` to `(nothing, nothing)`
        d = Dict()
        e_map = Dict()
        for net in constraints
            for e in net.edge
                if e.hybrid
                    if !e.ismajor
                        d[getchild(e)] = Vector{EdgeOrNA}([nothing, nothing, nothing])      # first entry:  `from` edge,
                                                                                            # second entry: `to` edge,
                                                                                            # third entry:  second `to` edge (i.e. if a node has 2 hybrid children)
                    end
                    e_map[e] = getchild(e)
                end
            end
        end
        new(d, e_map)
    end
end


function fully_logged(r::ReticMap, hyb::Node)
    return r.map[hyb][1] !== nothing && (r.map[hyb][2] !== nothing || r.map[hyb][3] !== nothing)
end


# Try logging a retic to "from" or "to", being lenient of errors (used exclusively for rootretic updating)
function trylogretic!(r::ReticMap, hyb::Node, subnetedge::Edge, fromorto::String)
    try
        logretic!(r, hyb, subnetedge, fromorto)
    catch e
        try
            logretic!(r, hyb, subnetedge, fromorto == "from" ? "to" : "from")
        catch
        end
    end
end
trylogretic!(r::ReticMap, hyb_edge::Edge, subnetedge::Edge, fromorto::String) =
    trylogretic!(r, r.edge_map[hyb_edge], subnetedge, fromorto)

function trylogretic_single!(r::ReticMap, hyb::Node, subnetedge::Edge, fromorto::String)
    try
        logretic!(r, hyb, subnetedge, fromorto)
    catch e
    end
end
trylogretic_single!(r::ReticMap, hyb_edge::Edge, subnetedge::Edge, fromorto::String) =
    trylogretic_single!(r, r.edge_map[hyb_edge], subnetedge, fromorto)

function logretic!(r::ReticMap, hyb::Node, subnetedge::Edge, fromorto::String)
    # If we're double logging identical edges then return w/o error
    if haskey(r.map, hyb)
        if fromorto == "from" && r.map[hyb][1] == subnetedge return end
        if fromorto == "to" && r.map[hyb][2] == subnetedge return end
    end

    # Log the reticulation
    if fromorto == "from"
        if (r.map[hyb][1] !== nothing)
            if r.map[hyb][1] != subnetedge
                throw(ErrorException("Overriding `from` edge"))
            end
        end
        if r.map[hyb][2] == subnetedge
            throw(ErrorException("Attempting to set `from` edge to a duplicate of the `to` edge."))
        elseif r.map[hyb][3] == subnetedge
            throw(ErrorException("Attempting to set `from` edge to a duplicate of the (second) `to` edge."))
        end
        r.map[hyb][1] = subnetedge
    else
        if (r.map[hyb][2] !== nothing)
            if (r.map[hyb][3] !== nothing)
                throw(ErrorException("Both `to` edges already set to non-nothing values."))
            else
                if r.map[hyb][1] == subnetedge
                    throw(ErrorException("Attempting to set (second) `to` edge to a duplicate of the `from` edge."))
                end
                r.map[hyb][3] = subnetedge
            end
        else
            if r.map[hyb][1] == subnetedge
                throw(ErrorException("Attempting to set `to` edge to a duplicate of the `from` edge."))
            end
            r.map[hyb][2] = subnetedge
        end
    end
end
logretic!(r::ReticMap, hyb_edge::Edge, subnetedge::Edge, fromorto::String) =
    logretic!(r, r.edge_map[hyb_edge], subnetedge, fromorto)

function check_reticmap(r::ReticMap)
    for (i, key) in enumerate(keys(r.map))
        if length(r.map[key]) != 3
            throw(ErrorException("ReticMap key $i has $(length(r.map[key])) attached edges."))
        elseif sum(r.map[key] .!== nothing) != 2 && sum(r.map[key] .!== nothing) != 3
            println("r.map length: $(length(r.map))")
            println(r.map[key])
            println(key.number)
            throw(ErrorException("ReticMap key $i has $(sum(r.map[key] .!== nothing)) attached non-nothing edges."))
        elseif r.map[key][1] === nothing
            @show key.name
            throw(ErrorException("ReitcMap key $i has no \"from\" edge."))
        end
    end
end
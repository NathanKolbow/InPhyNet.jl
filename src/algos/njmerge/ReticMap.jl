const EdgeOrNA = Union{Edge, Nothing}
struct ReticMap
    map::Dict{Edge, Tuple{EdgeOrNA, EdgeOrNA}}
    ReticMap() = error("not yet implemented")
end
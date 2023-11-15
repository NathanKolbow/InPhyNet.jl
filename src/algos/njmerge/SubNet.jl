using PhyloNetworks

struct SubNet
    # Fields
    nodes::Vector{Node}
    edges::Vector{Edge}
    linkpoint::Node
    id::Int64

    # Constructors
    SubNet() = error("not yet implemented")
    function SubNet(id, name::AbstractString)
        taxanode = Node(0, true)
        taxanode.name = name

        new([taxanode], [], taxanode, id)
    end
    SubNet(nodes::Vector{Node}, edges::Vector{Edge}, link::Node, newid::Int64) = new(nodes, edges, link, newid)
end


"""

Merges two unrooted SubNets into a single unrooted SubNet
Returns the merged subnet AND the edges in the merged net
    corresponding to where reticulations will be placed if there
    are corresponding reticulations discovered in the constraints.

IMPORTANT: this function should NOT make copies of edges or
    nodes, object references are stored and used in
    the main merging algorithm
"""
function mergesubnets!(n1::SubNet, n2::SubNet; newid::Int64=n1.id)
    newlink = Node(0, false)
    e1 = connectnodes!(newlink, n1.linkpoint)
    e2 = connectnodes!(newlink, n2.linkpoint)

    return SubNet([n1.nodes; n2.nodes; newlink], [n1.edges; n2.edges; e1; e2], newlink, newid), e1, e2
end


"""
Deals with all the overhead of connecting two nodes with an edge.
Returns the edge used to connect the nodes.

TODO: add option to make the connection a hybrid connection
"""
function connectnodes!(child::Node, parent::Node)
    edge = Edge(0, -1.)
    edge.node = [child, parent]
    push!(child.edge, edge)
    push!(parent.edge, edge)
    edge.containRoot = false
    edge.isChild1 = true
    return edge
end
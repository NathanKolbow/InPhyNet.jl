# Contains info for which species/taxa/etc. this lineage is an ancestor for
# e.g. if the bit array `lineage` is [0, 0, 0, 1, 1, 0], then this lineage
# is an ancestor of species 4 & 5, but not of any others
struct Lineage
    # Fields
    lineage::BitArray
    
    # Constructors
    function Lineage(i::Real, ntaxa::Real)
        Lineage([i], ntaxa)
    end
    function Lineage(members::AbstractVector{<:Real}, ntaxa::Real)
        lineage = BitArray([i in members for i=1:ntaxa])
        new(lineage)
    end
    function Lineage(l1::Lineage, l2::Lineage)
        # ntaxa(l1) == ntaxa(l2) || error("l1 and l2 have different lengths.")
        new(l1.lineage .|| l2.lineage)
    end
    function Lineage()
        new([])
    end
end

# Helper methods
ntaxa(lineage::Lineage) = length(lineage.lineage)
Base.coalesce(l1::Lineage, l2::Lineage) = Lineage(l1, l2)
Base.coalesce(ls::AbstractVector{Lineage}) = reduce(coalesce, ls)

const LDict = Dict{PhyloNetworks.Node, Union{LineageNode, Nothing}}
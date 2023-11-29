# Used to track tree probabilities, and all 
# LineageNodes associated with a given network
struct InterimParentalTree
    # Fields
    net::HybridNetwork  # does not necessarily reflect the actual topology of the IPT
    p::BigFloat         # probability *thus far* for this IPT

    # Constructors
    InterimParentalTree(net::HybridNetwork) = new(net, BigFloat(1.))
    InterimParentalTree(net::HybridNetwork, p::Real) = new(net, BigFloat(p))
end
const IPT = InterimParentalTree

# Helpers
top(t::IPT) = t.net
p(t::IPT) = t.p
prob(t::IPT) = t.p
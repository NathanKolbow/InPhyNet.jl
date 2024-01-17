"""
Calculates the `nretics_inside`, `nretics_outside`, and `nretics_duplicated` metrics.
"""
function calculateReticData(truenet::HybridNetwork, constraints::Vector{HybridNetwork})
    true_retic_names = [retic.name for retic in truenet.hybrid]
    constraint_retic_names = Set{String}()

    nretics_duplicated = 0
    for c in constraints
        # Log retics, counting duplicates
        c_reticnames = [retic.name for retic in c.hybrid]
        for retic in c_reticnames
            if retic in constraint_retic_names
                nretics_duplicated += 1
            else
                push!(constraint_retic_names, retic)
            end
        end
    end

    nretics_inside = length(constraint_retic_names)
    nretics_outside = truenet.numHybrids - nretics_inside

    return nretics_inside, nretics_outside, nretics_duplicated
end
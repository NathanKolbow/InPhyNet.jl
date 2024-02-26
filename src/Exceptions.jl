
struct SolutionDNEError <: Exception end
Base.showerror(io::IO, e::SolutionDNEError) = print(io, "No valid (i, j) pairs exist for the given distance matrix.")

struct ConstraintError <: Exception
    idx::Int64
    errmsg::String
end
Base.showerror(io::IO, e::ConstraintError) = print(io, "Constraint #$(e.idx) is invalid: $(e.errmsg)")
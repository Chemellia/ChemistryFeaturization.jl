# Type to store atomic graphs
struct AtomGraph{G <: AbstractGraph}
    graph::G
    elements::Vector{String}
    lapl::Matrix{Float32}
    features::Matrix{Float32}
    feat::Ref{Vector{AtomFeat}}
end

# define constructors, check if things are the right lengths

# define some functions
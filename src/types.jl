# Type to store featurization metadata. An array of these specifies a featurization scheme.
struct AtomFeat
    name::String
    categorial::Bool
    bins::Int
    logspaced::Bool
    vals::Array
end

# define constructors

# define some functions specialized on whether it's categorial or not, etc.


# Type to store atomic graphs
struct AtomGraph{G <: AbstractGraph}
    graph::G
    elements::Vector{String}
    lapl::Matrix{Float32}
    feature::Matrix{Float32}
end

# define constructors, check if things are the right lengths

# define some functions
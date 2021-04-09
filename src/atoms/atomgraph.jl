# some commentary

# docstring
mutable struct AtomGraph <: AbstractAtoms
    graph::SimpleWeightedGraph{Integer,Real}
    elements::Vector{String}
    lapl::LightGraphs.LinAlg.NormalizedLaplacian
    atom_feats::Matrix{Real} # if we add edge features this type will have to relax
    featurization::GraphNodeFeaturization
    id::String # or maybe we let it be a number too?
end

# constructors etc. as before
# maybe some cutesy stuff like dispatching things like length as length(elements)
# pretty printing stuff

# some commentary

struct GraphNodeFeaturization <: AbstractFeaturization
    atom_feats::Vector{AtomFeat}
    feature_vectors::Dict{String,Vector{Real}} # map from element symbol to vector
end

# docstring
function featurize!(a::AtomGraph, f::GraphNodeFeaturization)
    # if adding the exclude option to generic fcn, may not need to dispatch here
end
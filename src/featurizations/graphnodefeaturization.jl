#=
Featurization for `AtomGraph` objects that featurizes graph nodes with `AtomFeat` objects that depend only on elemental identity (i.e. `contextual==false`), and hence can be pre-tabulated.
=#
struct GraphNodeFeaturization <: AbstractFeaturization
    atom_feats::Vector{AtomFeature}
    feature_vectors::Dict{String,Vector{Real}} # map from element symbol to vector
    # add default constructor that checks that contextual==true, `feature_vectors` are of correct length, keys are elements, etc.
end

# TODO: in utils, add fcns to construct from feature names, etc.

# docstring
function featurize!(a::AtomGraph, f::GraphNodeFeaturization)
    # loop over nodes and just pull from `feature_vectors`
end
# TODO: proper docstring
#=
Featurization for `AtomGraph` objects that featurizes graph nodes only.
=#
struct GraphNodeFeaturization <: AbstractFeaturization
    atom_feats::Vector{AtomFeature}
end

# TODO: function to compute total vector length

# TODO: in utils, add fcns to construct from feature names, etc.
# (needs to build feature vectors also)
# TODO: proper docstring
#=
Featurization for `AtomGraph` objects that featurizes graph nodes only.
=#
struct GraphNodeFeaturization <: AbstractFeaturization
    atom_features::Vector{AtomFeature}
end

# TODO: pretty printing

# TODO: in utils, add fcns to construct from feature names, etc.
# (needs to build feature vectors also)
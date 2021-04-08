#= FEATURIZATION OBJECTS
All such objects should define at least one list of <:AbstractFeature objects and either work according to the generic featurize! defined herein or dispatch featurize! if customized behavior is needed.
=#

module Featurizations

# include...
include("graphnodefeaturization.jl")
include("weavefeaturization.jl")

# export...


# maybe also define abstract type here
abstract type AbstractFeaturization end

# generic featurize...there's probably a better way to write this...
# TODO: rewrite Dhairya's way
# docstring
function featurize!(a<:AbstractAtoms, f<:AbstractFeaturization)
    feats_lists = [field for field in fieldnames(typeof(f)) if field!==:combine]
    encoded_features = Dict(feats_list => [] for feats_list in feats_lists)
    for feats_list in feats_lists
        features = map(feat->feat(a), feats_list)
        encoded_features[feats_list] = features
    end
    a.features = f.combine(encoded_features)
    a.featurization = f
end

end

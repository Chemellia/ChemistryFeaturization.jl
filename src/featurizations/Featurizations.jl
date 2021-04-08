module Featurizations

# include...
include("graphnodefeaturization.jl")
include("weavefeaturization.jl")

# export...


# maybe also define abstract type here
abstract type AbstractFeaturization end

# generic featurize...there's probably a better way to write this...
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

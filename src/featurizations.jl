"""
    AbstractFeaturization

 A featurization stores a set of FeatureDescriptors and associated Codecs and defines how to combine the encoded values into whatever format is required to feed into a model.
"""
abstract type AbstractFeaturization end

"""
    features(featurization::AbstractFeaturization)

Return the list of feature descriptors used by `featurization`.
"""
features(fzn::AbstractFeaturization)::Vector{<:AbstractFeatureDescriptor} =
    throw(MethodError(features, fzn))


"""
    encodable_elements(featurization)

Return a list of elemental symbols that are valid constituents for structures that `featurization` can featurize.
"""
encodable_elements(fzn::AbstractFeaturization) =
    intersect(encodable_elements.(features(fzn))...)

"""
    encode(atoms, featurization)

Encode the features of `atoms` according to the scheme described by `featurization`.
"""
encode(atoms, fzn::AbstractFeaturization) = encode.(Ref(atoms), features(fzn))

"""
    decode(encoded, featurization)

Decode `encoded`, presuming it was encoded by `featurization`.
"""
decode(encoded, fzn::AbstractFeaturization) = decode.(encoded, features(fzn))

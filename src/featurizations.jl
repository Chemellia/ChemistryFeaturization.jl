"""
    AbstractFeaturization

 A featurization stores a set of FeatureDescriptors and associated Codecs and defines how to combine the encoded values into whatever format is required to feed into a model.
"""
abstract type AbstractFeaturization end

"""
    features(featurization)

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
    encode(featurization, atoms)

Encode the features of `atoms` according to the scheme described by `featurization`.
"""
encode(fzn::AbstractFeaturization, atoms) = encode.(features(fzn), Ref(atoms))

"""
    decode(featurization, encoded)

Decode `encoded`, presuming it was encoded by `featurization`.
"""
decode(fzn::AbstractFeaturization, encoded) = decode.(features(fzn), encoded)

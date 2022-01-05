module ChemistryFeaturization

# encodable_elements(a::Any) = throw(MethodError(encodable_elements, a))
# encode(a::Any, object_to_be_encoded) = throw(MethodError(encode, a))
# decode(a::Any, encoded_features) = throw(MethodError(decode, a))
# elements(a::Any) = throw(MethodError(elements, a))
# output_shape(a::Any) = throw(MethodError(output_shape, a))

include("data.jl")
export Data

include("codecs/codecs.jl")
export AbstractCodec, SimpleCodec, OneHotOneCold, DirectCodec
export encode, decode, output_shape

include("features/features.jl")
export AbstractFeatureDescriptor, 
    AbstractAtomFeatureDescriptor, 
    AbstractPairFeatureDescriptor
export get_value, default_codec

include("features/elementfeature.jl")
export ElementFeature

include("features/pairfeature.jl")


include("featurizations.jl")

include("featurizedatoms.jl")
export FeaturizedAtoms, featurize

end

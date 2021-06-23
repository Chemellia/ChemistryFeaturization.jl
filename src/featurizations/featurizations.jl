module Featurization

using ..ChemistryFeaturization.AbstractType: AbstractFeaturization, AbstractAtoms, AbstractCodec

import ..ChemistryFeaturization.encodable_elements
encodable_elements(fzn::AbstractFeaturization) = throw(MethodError(encodable_elements, fzn))
export encodable_elements

featurize!(a::AbstractAtoms, fzn::AbstractFeaturization) =
    throw(MethodError(featurize!, (a, fzn)))

import ..ChemistryFeaturization.decode
decode(fzn::AbstractFeaturization, encoded_feature) =
    throw(MethodError(decode, (fzn, encoded_feature)))
include("graphnodefeaturization.jl")
export GraphNodeFeaturization, featurize!, decode

validate_features(fzn::AbstractFeaturization, ed::AbstractCodec, encoded) = throw(MethodError(validate_features, (fzn, ed, encoded)))
export validate_features

include("weavefeaturization.jl")
export WeaveFeaturization

include("dummyfeaturization.jl")
export DummyFeaturization

end

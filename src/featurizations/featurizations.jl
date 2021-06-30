module Featurization

using ..ChemistryFeaturization.AbstractType: AbstractFeaturization

import ..ChemistryFeaturization.encodable_elements
encodable_elements(fzn::AbstractFeaturization) = throw(MethodError(encodable_elements, fzn))
export encodable_elements

import ..ChemistryFeaturization.encode
encode(fzn::AbstractFeaturization, object_to_be_encoded) = throw(MethodError(encode, fzn))
export encode

import ..ChemistryFeaturization.decode
decode(fzn::AbstractFeaturization, encoded_feature) = throw(MethodError(decode, fzn))
export decode

include("graphnodefeaturization.jl")
export GraphNodeFeaturization, encode, decode

include("weavefeaturization.jl")
export WeaveFeaturization

end

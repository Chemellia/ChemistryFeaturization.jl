module Featurization

using ..ChemistryFeaturization.AbstractType: AbstractFeaturization

import ..ChemistryFeaturization.encodable_elements
encodable_elements(fzn::AbstractFeaturization) = throw(MethodError(encodable_elements, fzn))
export encodable_elements

import ..ChemistryFeaturization.decode
decode(fzn::AbstractFeaturization, encoded_feature) = throw(MethodError(decode, fzn))
include("graphnodefeaturization.jl")
export GraphNodeFeaturization, featurize!, decode

include("weavefeaturization.jl")
export WeaveFeaturization

include("dummyfeaturization.jl")
export DummyFeaturization

end

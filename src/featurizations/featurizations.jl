module Featurization

using ..ChemistryFeaturization.AbstractType: AbstractFeaturization, AbstractAtoms

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

include("weavefeaturization.jl")
export WeaveFeaturization

include("dummyfeaturization.jl")
export DummyFeaturization

end

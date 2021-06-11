module Featurizations

using ..ChemistryFeaturization.AbstractTypes: AbstractFeaturization

import ..ChemistryFeaturization.encodable_elements
encodable_elements(fzn::AbstractFeaturization) = throw(MethodError(encodable_elements, fzn))
export encodable_elements

include("graphnodefeaturization.jl")
export GraphNodeFeaturization, featurize!, decode

include("weavefeaturization.jl")
export WeaveFeaturization

end

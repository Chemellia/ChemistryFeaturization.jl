module Featurizations

# basically export everything but `chunk_vec`
include("graphnodefeaturization.jl")
export GraphNodeFeaturization, encodable_elements, featurize!, decode

include("weavefeaturization.jl")
export WeaveFeaturization, encodable_elements

end

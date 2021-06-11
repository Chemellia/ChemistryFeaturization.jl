module ChemistryFeaturization

using SimpleWeightedGraphs
using Reexport

include("utils/Utils.jl")
export Utils
@reexport using .Utils.ElementFeatureUtils
@reexport using .Utils.GraphBuilding

include("abstracts/abstracttypes.jl")
export AbstractTypes
include("codecs/codecs.jl")
export Codecs

encodable_elements(a::Any) = throw(MethodError(encodable_elements, a))
decode(a::Any, encoded_features) = throw(MethodError(decode, a))

include("features/features.jl")
export FeatureDescriptors
include("atoms/atoms.jl")
export Atoms
include("featurizations/featurizations.jl")
export Featurizations

export encodable_elements, decode

end

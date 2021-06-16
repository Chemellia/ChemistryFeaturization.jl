module ChemistryFeaturization

using SimpleWeightedGraphs
using Reexport

include("utils/Utils.jl")
export Utils
@reexport using .Utils.ElementFeatureUtils
@reexport using .Utils.GraphBuilding

include("abstracts/abstracttypes.jl")
export AbstractType

include("codecs/codecs.jl")
export Codec

encodable_elements(a::Any) = throw(MethodError(encodable_elements, a))
decode(a::Any, encoded_features) = throw(MethodError(decode, a))

include("features/features.jl")
export FeatureDescriptor
export ElementFeatureDescriptor

include("atoms/atoms.jl")
export Atoms
export AtomGraph

include("featurizations/featurizations.jl")
export Featurization
export GraphNodeFeaturization

export encodable_elements, decode

end

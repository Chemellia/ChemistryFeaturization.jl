module ChemistryFeaturization

using SimpleWeightedGraphs

encodable_elements(a::Any) = throw(MethodError(encodable_elements, a))
encode(a::Any, object_to_be_encoded) = throw(MethodError(encode, a))
decode(a::Any, encoded_features) = throw(MethodError(decode, a))

include("utils/Utils.jl")
export Utils

include("abstracts/abstracttypes.jl")
export AbstractType

include("codecs/codecs.jl")
export Codec

include("features/features.jl")
export FeatureDescriptor
using .FeatureDescriptor: ElementFeatureDescriptor
export ElementFeatureDescriptor

include("atoms/atoms.jl")
export Atoms

using .Atoms: AtomGraph, visualize
export AtomGraph, visualize

include("featurizations/featurizations.jl")
export Featurization
using .Featurization: GraphNodeFeaturization, encode
export GraphNodeFeaturization, encode

export encodable_elements, decode

include("featurizedatoms.jl")
export FeaturizedAtoms

end

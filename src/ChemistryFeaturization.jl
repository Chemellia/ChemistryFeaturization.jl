module ChemistryFeaturization

using SimpleWeightedGraphs

encodable_elements(a::Any) = throw(MethodError(encodable_elements, a))
encode(a::Any, object_to_be_encoded) = throw(MethodError(encode, a))
decode(a::Any, encoded_features) = throw(MethodError(decode, a))
elements(a::Any) = throw(MethodError(elements, a))
output_shape(a::Any) = throw(MethodError(output_shape, a))

include("data.jl")
export Data

include("abstracts/abstracttypes.jl")
export AbstractType

include("codecs/codecs.jl")
export Codec

include("utils/Utils.jl")
export Utils

include("atoms/atoms.jl")
export Atoms

using .Atoms: AtomGraph, visualize, elements
export AtomGraph, visualize

include("features/features.jl")
export FeatureDescriptor
using .FeatureDescriptor
export ElementFeatureDescriptor, SpeciesFeatureDescriptor, BondFeatureDescriptor
export output_shape, get_value

include("featurizations/featurizations.jl")
export Featurization
using .Featurization: GraphNodeFeaturization, WeaveFeaturization, encode
export GraphNodeFeaturization, WeaveFeaturization, encode

export encodable_elements, decode

include("featurizedatoms.jl")
export FeaturizedAtoms, featurize, decode

end

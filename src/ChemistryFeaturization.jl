module ChemistryFeaturization
using AtomsBase

export elements
"""
    elements(atoms)

Return the list of elemental symbols corresponding to the atoms making up `atoms`.
"""
elements(atoms) = throw(MethodError(elements, atoms))
elements(sys::AbstractAtomicSystem) = String.(atomic_symbol(sys))

include("data.jl")
export Data

include("codecs/codecs.jl")
export AbstractCodec, SimpleCodec, OneHotOneCold, DirectCodec
export encode, decode, output_shape

include("features/features.jl")
export AbstractFeatureDescriptor, 
    AbstractAtomFeatureDescriptor, 
    AbstractPairFeatureDescriptor
export get_value, default_codec, encodable_elements

include("features/elementfeature.jl")
export ElementFeature

include("featurizations.jl")
export AbstractFeaturization, features

include("featurizedatoms.jl")
export FeaturizedAtoms, featurize

end

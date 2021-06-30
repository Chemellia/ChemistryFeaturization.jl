module FeatureDescriptor

using ..ChemistryFeaturization.AbstractType: AbstractFeatureDescriptor

import ..ChemistryFeaturization.encodable_elements
encodable_elements(fd::AbstractFeatureDescriptor) =
    throw(MethodError(encodable_elements, fd))
export encodable_elements

import ..ChemistryFeaturization.encode
encode(fd::AbstractFeatureDescriptor, object_to_be_encoded) = throw(MethodError(encode, fd))
export encode

import ..ChemistryFeaturization.decode
decode(fd::AbstractFeatureDescriptor, encoded_feature) = throw(MethodError(decode, fd))
export decode

include("abstractfeatures.jl")

include("bondfeatures.jl")

include("elementfeature.jl")
export ElementFeatureDescriptor, encode, decode, output_shape

include("ofm.jl")
export OrbitalFieldMatrix

include("pairfeature.jl")
export PairFeatureDescriptor

include("speciesfeature.jl")
export SpeciesFeatureDescriptor

end

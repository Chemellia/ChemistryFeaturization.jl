module FeatureDescriptor

using ..ChemistryFeaturization.AbstractType: AbstractFeatureDescriptor, AbstractCodec

import ..ChemistryFeaturization.encodable_elements
encodable_elements(fd::AbstractFeatureDescriptor) =
    throw(MethodError(encodable_elements, fd))
export encodable_elements

import ..ChemistryFeaturization.decode
decode(fd::AbstractFeatureDescriptor, encoded_feature) = throw(MethodError(decode, fd))
export decode

output_shape(fd::AbstractFeatureDescriptor, ed::AbstractCodec) = throw(MethodError(output_shape, (fd, ed)))
output_shape(fd::AbstractFeatureDescriptor) = output_shape(fd, fd.encoder_decoder)
export output_shape

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

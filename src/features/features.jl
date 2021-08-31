"""
    FeatureDescriptor

A feature, in this context, refers to a distinct attribute manifested
by an atomic system.
A FeatureDescriptor describes a feature/class of features - i.e., its name,
possible values, etc.

A FeatureDescriptor does NOT store any actual instances of the value(s) of the
feature it describes.
Simply put, it can be understood to be "features of a feature".

All FeatureDescriptors MUST also describe an encoding and decoding scheme.
This can (and should) be easily done using a Codec.
"""
module FeatureDescriptor

using Base: Int16
using ..ChemistryFeaturization.AbstractType:
    AbstractAtoms, AbstractFeatureDescriptor, AbstractCodec

import ..ChemistryFeaturization.encodable_elements
encodable_elements(fd::AbstractFeatureDescriptor) =
    throw(MethodError(encodable_elements, fd))
export encodable_elements

import ..ChemistryFeaturization.encode
"""
    encode(fd::AbstractAtomFeatureDescriptor, atoms::AbstractAtoms)
Encode `atoms` using the feature descriptor `fd`.
"""
encode(fd::AbstractFeatureDescriptor, atoms::AbstractAtoms) = fd(atoms)
export encode

import ..ChemistryFeaturization.decode
"""
    decode(fd::AbstractAtomFeatureDescriptor, encoded_feature)
Decode `encoded_feature` using the feature descriptor `fd`.
"""
decode(fd::AbstractFeatureDescriptor, encoded_feature) =
    fd.encoder_decoder(fd, encoded_feature)
export decode

(codec::AbstractCodec)(fd::AbstractFeatureDescriptor, a::AbstractAtoms) = error(
    "Logic specifying how $(typeof(codec))'s encoding mechanism actually encodes $(typeof(fd)) is undefined.",
)
(codec::AbstractCodec)(fd::AbstractFeatureDescriptor, encoded_feature) = error(
    "Logic specifying how $(typeof(codec))'s decoding mechanism actually decodes $(typeof(fd)) is undefined.",
)

output_shape(efd::AbstractFeatureDescriptor) = output_shape(fd, fd.encoder_decoder)
export output_shape

include("abstractfeatures.jl")

include("bondfeatures.jl")

include("elementfeature.jl")
export ElementFeatureDescriptor, encode, decode, output_shape

include("orbitalfeature.jl")
export OrbitalFeatureDescriptor

include("pairfeature.jl")
export PairFeatureDescriptor

include("speciesfeature.jl")
export SpeciesFeatureDescriptor

end

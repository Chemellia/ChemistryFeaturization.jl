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

"""
    get_value(fd::AbstractFeatureDescriptor, atoms::AbstractAtoms)
Get the value(s) of feature corresponding to feature descriptor `fd` for structure `atoms`.

See also: [`encode`](@ref)
"""
get_value(fd::AbstractFeatureDescriptor, atoms::AbstractAtoms) = throw(MethodError(fd, atoms))
(fd::AbstractFeatureDescriptor)(atoms::AbstractAtoms) = get_value(fd, atoms)

import ..ChemistryFeaturization.encode
"""
    encode(fd::AbstractFeatureDescriptor, atoms::AbstractAtoms)
Encode features for `atoms` using the feature descriptor `fd`.
"""
encode(fd::AbstractFeatureDescriptor, atoms::AbstractAtoms) = encode(fd.encoder_decoder, get_value(fd, atoms)) # TODO: would like this to be get_value for clarity but scope issues at present
export encode

import ..ChemistryFeaturization.decode
"""
    decode(fd::AbstractFeatureDescriptor, encoded_feature)
Decode `encoded_feature` using the feature descriptor `fd`.
"""
decode(fd::AbstractFeatureDescriptor, encoded_feature) = decode(fd.encoder_decoder, encoded_feature)
export decode

output_shape(efd::AbstractFeatureDescriptor) = output_shape(efd, efd.encoder_decoder)
export output_shape

include("abstractfeatures.jl")
encode(efd::AbstractAtomFeatureDescriptor, atoms::AbstractAtoms) = hcat(encode(efd.encoder_decoder, get_value(efd, atoms))...)

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

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

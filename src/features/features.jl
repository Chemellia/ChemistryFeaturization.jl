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

using Base: Int16
using .Codec: OneHotOneCold

"""
    AbstractFeatureDescriptor

All [FeatureDescriptors](@ref fd) defined for different types of features must be a
subtype of AbstractFeatureDescriptor.
"""
abstract type AbstractFeatureDescriptor end

# pretty printing, short version
Base.show(io::IO, fd::AbstractFeatureDescriptor) = print(io, "$(typeof(fd)) $(fd.name)")

encodable_elements(fd::AbstractFeatureDescriptor) =
    throw(MethodError(encodable_elements, fd))

"""
    get_value(fd::AbstractFeatureDescriptor, atoms)
Get the value(s) of feature corresponding to feature descriptor `fd` for structure `atoms`.
This function computes and returns the value that would actually get encoded by [`encode`](@ref).
"""
get_value(fd::AbstractFeatureDescriptor, atoms) =
    throw(MethodError(fd, atoms))
(fd::AbstractFeatureDescriptor)(atoms) = get_value(fd, atoms)

"""
    encode(fd::AbstractFeatureDescriptor, atoms)
Encode features for `atoms` using the feature descriptor `fd`.
"""
encode(fd::AbstractFeatureDescriptor, atoms) =
    encode(fd.encoder_decoder, get_value(fd, atoms))
export encode
export get_value

"""
    decode(fd::AbstractFeatureDescriptor, encoded_feature)
Decode `encoded_feature` using the feature descriptor `fd`.
"""
decode(fd::AbstractFeatureDescriptor, encoded_feature) =
    decode(fd.encoder_decoder, encoded_feature)
export decode

export AbstractAtomFeatureDescriptor,
    AbstractPairFeatureDescriptor, AbstractEnvironmentFeatureDescriptor

abstract type AbstractAtomFeatureDescriptor <: AbstractFeatureDescriptor end
abstract type AbstractPairFeatureDescriptor <: AbstractFeatureDescriptor end

encode(efd::AbstractAtomFeatureDescriptor, atoms) =
    hcat(encode(efd.encoder_decoder, get_value(efd, atoms))...)

output_shape(afd::AbstractFeatureDescriptor) = output_shape(afd, afd.encoder_decoder)
output_shape(::AbstractAtomFeatureDescriptor, ed::OneHotOneCold) = output_shape(ed)

include("elementfeature.jl")


include("pairfeature.jl")
export PairFeatureDescriptor

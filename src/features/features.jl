"""
    FeatureDescriptor

A feature, in this context, refers to a distinct attribute manifested
by an atomic system.
A FeatureDescriptor describes a feature/class of features - i.e., its name,
possible values, etc.

A FeatureDescriptor does NOT store any actual instances of the value(s) of the
feature it describes.
Simply put, it can be understood to be "features of a feature".
"""

"""
    AbstractFeatureDescriptor

All [FeatureDescriptors](@ref fd) defined for different types of features must be a
subtype of AbstractFeatureDescriptor.
"""
abstract type AbstractFeatureDescriptor end

# pretty printing, short version
Base.show(io::IO, fd::AbstractFeatureDescriptor) = print(io, "$(typeof(fd)) $(fd.name)")

"""
    default_codec(fd::AbstractFeatureDescriptor)

Return the default codec to use for encoding values of `fd`.
"""
default_codec(fd::AbstractFeatureDescriptor) = throw(MethodError(default_codec, fd))

"""
    encodable_elements(fd::AbstractFeatureDescriptor)
Return a list of elemental symbols for which the feature associated with `fd` is defined.
"""
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
    encode(fd, atoms)
    encode(fd, codec, atoms)
Encode features for `atoms` using the feature descriptor `fd` using the default codec for `fd`. if `codec` is not specified.
"""
encode(fd::AbstractFeatureDescriptor, atoms) =
    encode(default_codec(fd), get_value(fd, atoms))
encode(fd::AbstractFeatureDescriptor, codec::AbstractCodec, atoms) =
    encode(codec, get_value(fd, atoms))

"""
    decode(fd, encoded_feature)
    decode(fd, codec, encoded_feature)
Decode `encoded_feature` using the feature descriptor `fd`, presuming it was encoded via `fd`'s default codec if `codec` is not specified.
"""
decode(fd::AbstractFeatureDescriptor, encoded_feature) =
    decode(default_codec(fd), encoded_feature)

abstract type AbstractAtomFeatureDescriptor <: AbstractFeatureDescriptor end
abstract type AbstractPairFeatureDescriptor <: AbstractFeatureDescriptor end

# TODO: check that these work properly
encode(afd::AbstractAtomFeatureDescriptor, atoms) =
    hcat(encode.(Ref(default_codec(afd)), get_value(afd, atoms))...)
encode(afd::AbstractAtomFeatureDescriptor, codec::AbstractCodec, atoms) = 
    hcat(encode.(Ref(codec), get_value(afd, atoms))...)

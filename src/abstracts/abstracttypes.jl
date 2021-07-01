module AbstractType

"""
    AbstractAtoms

Atoms objects characterize structural representations of molecular structures.
These representations can either be standardized forms
(like graph representation using an adjacency matrix, point cloud representation, etc),
or be custom user-defined representations too.

All types defined for different representations of atomic systems
must be a subtype of AbstractAtoms.
"""
abstract type AbstractAtoms end


"""
    AbstractFeatureDescriptor

A feature, in this context, refers to a distinct attribute manifested
by an atomic system.
A FeatureDescriptor describes a feature/class of features - i.e., its name,
possible values, etc.

A FeatureDescriptor does NOT store any actual instances of the value(s) of the
feature it describes.
Simply put, it can be understood to be "features of a feature".

All FeatureDescriptors MUST also describe an encoding and decoding scheme.
This can (and should) be easily done using a Codec.

All FeatureDescriptors defined for different types of features must be a
subtype of AbstractFeatureDescriptor.
"""
abstract type AbstractFeatureDescriptor end


"""
    AbstractFeaturization

A Featurization collectively stores a set of FeatureDescriptors (and any other
supporting attributes that it may further require), and defines how to featurize
an Atoms object by using the encoding schemes defined by the FeatureDescriptors
stored.

All types defined for different featurizations must be a subtype
of AbstractFeaturization.
"""
abstract type AbstractFeaturization end

"""
    AbstractCodec

Codecs provide the encoding-decoding scheme in a flexible way.
Codecs are titled according to the general attributes that the codec scheme
would require, and NOT the actual encoding/decoding function itself.

Every Codec object MUST have one set of `encode_f::Function` and
`decode_f::Function` fields.

This makes things customizable, and allows plug-and-play behaviour with
different variants of the same codec scheme, or for strikingly similar codec
schemes.

All codecs defined for different encoding-decoding schemes
must be a subtype of AbstractCodec.
"""
abstract type AbstractCodec end

end

module AbstractType

"""
    AbstractAtoms

All types defined for different representations of atomic systems
(different [Atoms](@ref) types) must be a subtype of AbstractAtoms.
"""
abstract type AbstractAtoms end


"""
    AbstractFeatureDescriptor

All FeatureDescriptors defined for different types of features must be a
subtype of AbstractFeatureDescriptor.
"""
abstract type AbstractFeatureDescriptor end


"""
    AbstractFeaturization

All types defined for different featurization schemes must be a subtype
of AbstractFeaturization.
"""
abstract type AbstractFeaturization end

"""
    AbstractCodec

All codecs defined for different encoding-decoding schemes
must be a subtype of AbstractCodec.
"""
abstract type AbstractCodec end

end

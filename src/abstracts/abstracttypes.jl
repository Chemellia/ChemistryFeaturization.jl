module AbstractType

"""
    AbstractAtoms{A}

All types defined for different representations of atomic systems
(different [Atoms](@ref atoms) types) must be a subtype of AbstractAtoms. Type parameter represents the original structure representation from which the object was built.
"""
abstract type AbstractAtoms{A} end


"""
    AbstractFeatureDescriptor

All [FeatureDescriptors](@ref fd) defined for different types of features must be a
subtype of AbstractFeatureDescriptor.
"""
abstract type AbstractFeatureDescriptor end


"""
    AbstractFeaturization

All types defined for different [Featurizations](@ref fzn) must be a subtype
of AbstractFeaturization.
"""
abstract type AbstractFeaturization end

"""
    AbstractCodec

All [Codecs](@ref codecs) defined for different encoding-decoding schemes
must be a subtype of AbstractCodec.
"""
abstract type AbstractCodec end

end

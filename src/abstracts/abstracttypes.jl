module AbstractType

"""
    AbstractAtoms

Abstract type from which all concrete Atoms types are subtyped.
"""
abstract type AbstractAtoms end

"""
    AbstractFeatureDescriptor

Abstract type from which all concrete feature descriptor types are subtyped.
"""
abstract type AbstractFeatureDescriptor end
abstract type AbstractFeaturization end
abstract type AbstractCodec end

end

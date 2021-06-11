
using ..ChemistryFeaturization.AbstractTypes: AbstractFeatureDescriptor
export AbstractAtomFeatureDescriptor,
    AbstractPairFeatureDescriptor, AbstractEnvironmentFeatureDescriptor

abstract type AbstractAtomFeatureDescriptor <: AbstractFeatureDescriptor end
abstract type AbstractPairFeatureDescriptor <: AbstractFeatureDescriptor end
abstract type AbstractEnvironmentFeatureDescriptor <: AbstractFeatureDescriptor end

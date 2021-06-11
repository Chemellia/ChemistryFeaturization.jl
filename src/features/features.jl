module FeatureDescriptors

using ..ChemistryFeaturization.AbstractTypes: AbstractFeatureDescriptor

include("abstractfeatures.jl")

include("bondfeatures.jl")

include("elementfeature.jl")
# should these exports be put in the respective files itself?
export ElementFeatureDescriptor, encodable_elements, encode, decode

include("ofm.jl")
export OrbitalFieldMatrix

include("pairfeature.jl")
export PairFeatureDescriptor

include("speciesfeature.jl")
export SpeciesFeatureDescriptor, encodable_elements # will this work? (since encodable_elements is already being exported above...)

end

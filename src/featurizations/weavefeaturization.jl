# some commentary

struct WeaveFeaturization <: AbstractFeaturization
    atom_features::Vector{<:AbstractAtomFeatureDescriptor}
    pair_features::Vector{<:AbstractPairFeatureDescriptor}
end

function encodable_elements(fzn::WeaveFeaturization)
    # TODO: implement me!
end

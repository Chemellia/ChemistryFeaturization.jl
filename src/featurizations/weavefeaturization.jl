# some commentary

struct WeaveFeaturization <: AbstractFeaturization
    atom_features::Vector{<: AbstractAtomFeature}
    pair_features::Vector{<: AbstractPairFeature}
end

function encodable_elements(fzn::WeaveFeaturization)
    # TODO: implement me!
end
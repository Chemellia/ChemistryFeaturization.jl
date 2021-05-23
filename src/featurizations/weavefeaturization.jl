# some commentary

struct WeaveFeaturization <: AbstractFeaturization
    atom_features::Vector{AtomFeature}
    pair_features::Vector{GeneralPairFeature}
end

function encodable_elements(fzn::WeaveFeaturization)
    # TODO: implement me!
end
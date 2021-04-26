# some commentary

struct WeaveFeaturization <: AbstractFeaturization
    atom_features::Vector{AtomFeature}
    pair_features::Vector{PairFeature}
end
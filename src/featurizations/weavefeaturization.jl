# some commentary

struct WeaveFeaturization <: AbstractFeaturization
    atom_feats::Vector{AtomFeature}
    pair_feats::Vector{PairFeature}
end
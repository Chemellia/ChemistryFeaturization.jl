# some commentary

struct WeaveFeaturization <: AbstractFeaturization
    atom_feats::Vector{AtomFeat}
    pair_feats::Vector{PairFeat}
end
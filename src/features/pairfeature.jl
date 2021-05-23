#=
Composition of a GeneralPairFeature and BondFeature
=#
include("generalpairfeature.jl")
include("bondfeatures.jl")

# Maybe better as struct PairFeature{T, S} where T <: GeneralPairFeature, S <: BondFeature ?
struct PairFeature
    general_pair_feature::GeneralPairFeature
    bond_features::Union{Vector{BondFeature},Nothing}
end

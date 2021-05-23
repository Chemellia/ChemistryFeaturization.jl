#=
Composition of a PairFeature and BondFeature
=#
include("pairfeature.jl")
include("bondfeatures.jl")

# Maybe better as struct WeavePair{T, S} where T <: PairFeature, S <: BondFeature ?
struct WeavePair
    pair_feature::PairFeature
    bond_features::Union{Vector{BondFeature}, Nothing}
end

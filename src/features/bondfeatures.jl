include("abstractfeatures.jl")

abstract type BondFeatureDescriptor <: AbstractPairFeatureDescriptor end

#=
TODO - Come up with a better name for these type
        (since bond types can also be confused to be ionic, covalent, coordination bonds, etc.)
=#
# boolean features. existence of an object implies existence of feature.
abstract type Bonds end
struct Single <: Bonds end
struct Double <: Bonds end
struct Triple <: Bonds end

struct BondType{T<:Bonds} <: BondFeatureDescriptor end

struct InRing <: BondFeatureDescriptor end
struct IsConjugated <: BondFeatureDescriptor end

module ChemistryFeaturization

using SimpleWeightedGraphs
using Reexport

# define all the abstract types
export AbstractAtoms, AbstractFeatureDescriptor, AbstractFeaturization
export AbstractPairFeatureDescriptor,
    AbstractAtomFeatureDescriptor, AbstractEnvironmentFeatureDescriptor

abstract type AbstractAtoms end
abstract type AbstractFeatureDescriptor end
abstract type AbstractFeaturization end

abstract type AbstractAtomFeatureDescriptor <: AbstractFeatureDescriptor end
abstract type AbstractPairFeatureDescriptor <: AbstractFeatureDescriptor end
abstract type AbstractEnvironmentFeatureDescriptor <: AbstractFeatureDescriptor end

include("utils/Utils.jl")
@reexport using .Utils.ElementFeatureUtils
@reexport using .Utils.GraphBuilding

#= ATOMS OBJECTS

Should have a field called `elements` that is a list of strings of elemental symbols, as well as fields to store features and an associated featurization object.

AtomGraph shouldn't have to change too much from its current incarnation. WeaveMol is 
not currently a thing that can be directly fed into a model, so that will have to get
updated. Both will need some way to specify which types of features can be attached 
to them, maybe this is a place for the holy trait design pattern?

=#

# link to guidance in docs about how to implement new feature types

# export...
export AtomGraph, visualize
#export WeaveMol#, ...

# include...
include("atoms/atomgraph.jl")
#include("weavemol.jl")


#= FEATURE OBJECTS

I decided on an abstract type to allow the dispatch shown below in the "magical encoding"
and "magical decoding" bits.

Representative examples – feature: feature object type, input to encoder, output of encoder:

block (s,p,d,f):            ElementFeatureDescriptor, String, Flux.OneHotMatrix
electronegativity (binned): ElementFeatureDescriptor, Float32, Flux.OneHotMatrix
electronegativity (direct): ElementFeatureDescriptor, Float32, Vector{Float32}
oxidation state:            SpeciesFeatureDescriptor, Int, Vector{Float32}
distance between atoms:     PairFeat, Float32, Matrix{Float32}
bond type:                  PairFeat, String, Array{Float32,3}

(bond type is categorical and gets one-hot encoded, so the first two indices of the Matrix
should be atom indices and the third should be indexing into the one-hot vector, which should just contain all zeros if the two atoms are not bonded (or we add a bin to the one-hot encoding to indicate that)

All subtypes should have `encode_f` and `decode_f` fields
=#

# link to guidance in docs about how to implement new feature types

# export...
export ElementFeatureDescriptor, SpeciesFeatureDescriptor, PairFeatureDescriptor
export encodable_elements, decode

"""
    encodable_elements(f::AbstractFeatureDescriptor)
    encodable_elements(feature_name::String)
    encodable_elements(fzn::AbstractFeaturization)

Return a list of elements encodable by a given feature or featurization.
"""
encodable_elements(f::AbstractFeatureDescriptor) = println("Implement me please!")
encodable_elements(fzn::AbstractFeaturization) = println("Implement me please!")

# include...
include("features/elementfeature.jl")
include("features/speciesfeature.jl")
include("features/pairfeature.jl")

# generic encode
function (f::AbstractFeatureDescriptor)(a::AbstractAtoms)
    println("Implement me, please!")
end

"""
    decode(f::AbstractFeatureDescriptor, encoded_feature)

Decode a feature that was encoded using the provided feature object.

## Examples

```jldoctest
julia> decode(AtomFeatureDescriptor("Block"), [0, 1, 0, 0])
"p"
```
"""
function decode(f::AbstractFeatureDescriptor, encoded_feature)
    #f.decode_f(encoded_feature)
    println("Implement me, please!")
end

#= FEATURIZATION OBJECTS
All such objects should define at least one list of <: AbstractFeatureDescriptor objects and either work according to the generic featurize! defined herein or dispatch featurize! if customized behavior is needed. You should also dispatch the decode function, for which a generic implementation does not currently exist.
=#

# export...
export GraphNodeFeaturization, decode
export encodable_elements, featurize!
export WeaveFeaturization

# include...
include("featurizations/graphnodefeaturization.jl")
include("featurizations/weavefeaturization.jl")

"""
    featurize!(a::AbstractAtoms, fzn::AbstractFeaturization)

Featurize a structure with a given featurization.
"""
function featurize!(a::AbstractAtoms, fzn::AbstractFeaturization)
    println("Implement me please!")
end

# dispatch some things so that featurize!.(::Vector{AbstractAtoms}, ::AbstractFeaturization) will work...
Base.iterate(fzn::AbstractFeaturization) = (fzn, nothing)
Base.iterate(fzn::AbstractFeaturization, state) = nothing
Base.length(fzn::AbstractFeaturization) = 1

"""
    decode(fzn::AbstractFeaturization, encoded)

Decode a set of features that were encoded using the provided featurization.
"""
function decode(fzn::AbstractFeaturization, encoded)
    println("Implement me please!")
end

# pretty printing for featurizations, short version
Base.show(io::IO, fzn::AbstractFeaturization) = print(io, "$(typeof(fzn))")

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", fzn::AbstractFeaturization)
    st = "$(typeof(fzn)):"
    for feature_list in fieldnames(typeof(fzn))
        st = string(st, "\n    $(feature_list):")
        for feature in getproperty(fzn, feature_list)
            st = string(st, "\n        $(feature)")
        end
    end
    print(io, st)
end

end

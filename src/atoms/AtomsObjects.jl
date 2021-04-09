#= ATOMS OBJECTS

AtomGraph shouldn't have to change too much from its current incarnation. WeaveMol is 
not currently a thing that can be directly fed into a model, so that will have to get
updated. Both will need some way to specify which types of features can be attached 
to them, maybe this is a place for the holy trait design pattern?

Example: AtomGraph can take ElementFeat and ComputedAtomFeat but not PairFeat or StructureFeat, WeaveMol can take ElementFeat, ComputedAtomFeat, and PairFeat but also not StructureFeat
=#

module AtomsObjects

# link to guidance in docs about how to implement new feature types

# include...
include("atomgraph.jl")
include("weavemol.jl")

# export...
export AbstractAtoms, AtomGraph, WeaveMol

# maybe also define abstract type here
abstract type AbstractAtoms end

end

# [Feature Descriptors](@id fd)

```@docs
FeatureDescriptor
```

Feature descriptors store all necessary information to encode and decode feature values on various parts of an [Atoms](@ref atoms) object and appropriately combine them into a single form (vector, matrix, etc.), describing the value(s) of the feature for the entire object.

For example, if an `ElementFeatureDescriptor` encodes a vector for each atom in an object, they could be concatenated together into a matrix with a column for each atom to describe a structure.

Feature Descriptors must be designed with interoperability in mind. A [FeatureDescriptor](@ref fd) object must work deterministically with different datasets.

## Hierarchy

The type hierarchy of these objects is currently as follows.

```text
AbstractType.AbstractFeatureDescriptor
├─── AbstractAtomFeatureDescriptor
│    ├─── SpeciesFeatureDescriptor
│    └─── ElementFeatureDescriptor
│
├─── AbstractEnvironmentFeatureDescriptor
│    └─── OrbitalFieldMatrix
│
└─── AbstractPairFeatureDescriptor
     ├─── BondFeatureDescriptor
     │    ├─── BondType
     │    ├─── InRing
     │    └─── IsConjugated
     │
     └─── PairFeatureDescriptor
```

More details on each of these types can be found below.\
More types (e.g. environment features) will be implemented in the future!

## Functionality common to all feature descriptors

* they should be callable on atoms objects and return encoded features
* Similarly, `decode` should work...
```@docs
decode(::FeatureDescriptor.AbstractFeatureDescriptor, ::Any)
```
* the function `encodable_elements` should be defined on all feature descriptors, as it will be used to verify that a feature can be encoded for every atom in a structure.

## Atom Feature Descriptors

These types encode features for single atoms in a structure. The abstract parent type is `AtomFeatureDescriptor`.

### Element Feature Descriptors

An `ElementFeatureDescriptor`'s encoded values are defined only by the elemental identity of an atom. Examples include atomic mass and block (s, p, d, or f) in the periodic table.

```@docs
ElementFeatureDescriptor
```

In the example, below, we first show and then encode the value of the block of each atom in a hydrogen molecule. The result is two `hcat`ted vectors [1 0 0 0], indicating hydrogen is _s_-block.

```jldoctest; setup = :(using ChemistryFeaturization.Atoms, ChemistryFeaturization.FeatureDescriptor)
H2 = AtomGraph([0. 1.; 1. 0.], ["H", "H"])
block = ElementFeatureDescriptor("Block")
block(H2) # equivalent to get_value(block, H2)

# output
2-element Vector{String}:
 "s"
 "s"
```
```julia
julia> encode(block, H2)
4×2 Matrix{Float64}:
 1.0  1.0
 0.0  0.0
 0.0  0.0
 0.0  0.0
```

Because they are defined only by the element, values for these features can be tabulated in a lookup table. Many commonly-desired element features are included in the `atom_data_df` DataFrame, but you can also define custom lookup tables for other features by utilizing the `lookup_table` keyword of the `ElementFeatureDescriptor` constructor.

### Species Feature Descriptor

A `SpeciesFeatureDescriptor`'s encoded values depend on its local environment. Examples are an atom's formal oxidation state, or whether it is part of an aromatic ring.

```@docs
FeatureDescriptor.SpeciesFeatureDescriptor
```

A variety of built-in options for computing features on `GraphMol` objects (see MolecularGraph.jl package) are included. These can be constructed from a string, like

```julia
julia> sfd = SpeciesFeatureDescriptor("isaromatic")

SpeciesFeature isaromatic:
   categorical: true
   works on: GraphMol
```
This feature will return a an array of bits representing whether each atom in a structure is part of an aromatic ring. For example...
```julia
using MolecularGraph
using ChemistryFeaturization.Codec

caffeine = smilestomol("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
ag = AtomGraph(caffeine)
```
```julia
julia> get_value(sfd, ag)
14-element BitVector:
 0
 1
 1
 1
 1
 1
 1
 0
 1
 1
 0
 1
 0
 0
 ```

 We can also build custom `SpeciesFeatureDescriptors` "from scratch." For a simple example, let's suppose we want to encode graph-related features:
```julia
julia> ag = AtomGraph(Float32.([0 1; 1 1]), ["H", "O"]) # build a simple AtomGraph with a self loop
```
```julia
# a function that will take in a graph and return the number of neighbors of each node
num_nbs = g -> first.(length.(neighbors.(Ref(g), 1:nv(g))))
# we'll encode neighbors counts categorically in four bins
codec = OneHotOneCold(true, [1,2,3,4])
categorical = true
ee = ["O","H"]
# the type parameter should match the type parameter of the AtomGraph object
sfd = SpeciesFeatureDescriptor{SimpleWeightedGraph}("num_neighbors", num_nbs, codec, categorical, ee)
```

```julia
julia> get_value(sfd, ag)
2-element Vector{Int64}:
 1
 2
 ```

## Pair Feature Descriptors

Pair feature descriptors encode features of pairs of atoms. The abstract parent type is `AbstractPairFeatureDescriptor`.

The concrete type `PairFeatureDescriptor` encodes information about any pair of atoms, such as the distance between them.

### Bond Feature Descriptors

Bond feature descriptors are defined only for two atoms that are bonded to each other.

TODO: more details here
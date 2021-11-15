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

Bond feature descriptors are defined only for two atoms that are bonded to each other. Currently, the API is only built out for `AbstractAtoms` objects that have been constructed from `GraphMol` objects (defined by the MolecularGraph package). Similar to `SpeciesFeatureDescriptor`s, there are convenience functions defined in an associated utils file, and we can easily construct one like so:

```julia
julia> bfd = BondFeatureDescriptor("bondorder")
BondFeatureDescriptor{GraphMol, OneHotOneCold} bondorder:
   categorical: true
   works on: GraphMol
```

When we retrieve the values, note that they are now formatted as a (symmetric) matrix, with the two indices referring to the indices of the two atoms. If the atoms aren't bonded, that entry will be `missing`.
```julia
julia> ag = AtomGraph(smilestomol("c1ccncc1")) # pyridine
AtomGraph{GraphMol{SmilesAtom, SmilesBond}}  with 6 nodes, 6 edges
        atoms: ["C", "C", "C", "N", "C", "C"]
julia> get_value(bfd, ag)
6×6 Matrix{Union{Missing, Int64}}:
  missing  2          missing   missing   missing  1
 2          missing  1          missing   missing   missing
  missing  1          missing  2          missing   missing
  missing   missing  2          missing  1          missing
  missing   missing   missing  1          missing  2
 1          missing   missing   missing  2          missing
```

And when we encode it, the default result is now a 3D matrix (one vector for each pair of atoms). If the value is `missing`, the third dimension will be filled with `missing`:
```julia
julia> Array{Union{Missing,Int}}(encode(bfd, ag)) # cast to Int so printing looks nicer, you could leave this off but would just see trues and falses instead
6×6×3 Array{Union{Missing, Int64}, 3}:
[:, :, 1] =
  missing  0          missing   missing   missing  1
 0          missing  1          missing   missing   missing
  missing  1          missing  0          missing   missing
  missing   missing  0          missing  1          missing
  missing   missing   missing  1          missing  0
 1          missing   missing   missing  0          missing

[:, :, 2] =
  missing  1          missing   missing   missing  0
 1          missing  0          missing   missing   missing
  missing  0          missing  1          missing   missing
  missing   missing  1          missing  0          missing
  missing   missing   missing  0          missing  1
 0          missing   missing   missing  1          missing

[:, :, 3] =
  missing  0          missing   missing   missing  0
 0          missing  0          missing   missing   missing
  missing  0          missing  0          missing   missing
  missing   missing  0          missing  0          missing
  missing   missing   missing  0          missing  0
 0          missing   missing   missing  0          missing
```

Some of the built-in bond features are boolean valued, and by default are direct encoded (since onehot doesn't really make sense when there are only two options). In this case, the encoded matrix is equal to the values:

```julia
julia> bfd = BondFeatureDescriptor("isringbond")
BondFeatureDescriptor{GraphMol, DirectCodec} isringbond:
   categorical: true
   works on: GraphMol
julia> get_value(bfd, ag)
6×6 Matrix{Union{Missing, Bool}}:
     missing  true             missing      missing      missing  true
 true             missing  true             missing      missing      missing
     missing  true             missing  true             missing      missing
     missing      missing  true             missing  true             missing
     missing      missing      missing  true             missing  true
 true             missing      missing      missing  true             missing
 julia> encode(bfd, ag)
6×6 Matrix{Union{Missing, Int64}}:
  missing  1          missing   missing   missing  1
 1          missing  1          missing   missing   missing
  missing  1          missing  1          missing   missing
  missing   missing  1          missing  1          missing
  missing   missing   missing  1          missing  1
 1          missing   missing   missing  1          missing
```
(pyridine is a ring molecule, so every bond is a ring bond)
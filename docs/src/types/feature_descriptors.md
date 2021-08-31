# Feature Descriptors

Feature descriptors store all necessary information to encode **and decode** feature values on various parts of an atoms object and appropriately combine them into a single object (vector, matrix, etc.) describing the value/values of the feature for the entire object.

For example, if an `ElementFeatureDescriptor` encodes a vector for each atom in an object, they could be concatenated together into a matrix with a column for each atom to describe a structure.

The type hierarchy of these objects is currently:
```
|---- AbstractFeatureDescriptor
    |---- AbstractAtomFeatureDescriptor
        |==== ElementFeatureDescriptor
        |==== SpeciesFeatureDescriptor
    |---- AbstractPairFeatureDescriptor
        |==== PairFeatureDescriptor
        |---- BondFeatureDescriptor
            |==== BondType
            |==== InRing
            |==== IsConjugated
```
where 
`----` = Abstract Type
and
`====` = Concrete Type

More details on each of these types is below, and more types (e.g. environment features) will be implemented in the future!

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

In the example below, we encode the block of each atom in a hydrogen molecule. The result is two `hcat`-ed vectors [1 0 0 0], indicating hydrogen is _s_-block.

TODO: check that this test passes once new version is tagged

```jldoctest; setup = :(using ChemistryFeaturization.Atoms, ChemistryFeaturization.FeatureDescriptor)
H2 = AtomGraph([0. 1.; 1. 0.], ["H", "H"])
block = ElementFeatureDescriptor("Block")
block(H2)

# output
4×2 Matrix{Float64}:
 1.0  1.0
 0.0  0.0
 0.0  0.0
 0.0  0.0
```

Because they are defined only by the element, values for these features can be tabulated in a lookup table. Many commonly-desired element features are included in the `atom_data_df` DataFrame, but you can also define custom lookup tables for other features by utilizing the `lookup_table` keyword of the `ElementFeatureDescriptor` constructor.

TODO: add remark about encoding options once that PR is merged

### Species Feature Descriptor

A `SpeciesFeatureDescriptor`'s encoded values depend on its local environment. Examples are an atom's format oxidation state, or whether it is part of an aromatic ring.

TODO: more details once we have better examples

## Pair Feature Descriptors

Pair feature descriptors encode features of pairs of atoms. The abstract parent type is `AbstractPairFeatureDescriptor`.

The concrete type `PairFeatureDescriptor` encodes information about any pair of atoms, such as the distance between them.

### Bond Feature Descriptors

Bond feature descriptors are defined only for two atoms that are bonded to each other.

TODO: more details here
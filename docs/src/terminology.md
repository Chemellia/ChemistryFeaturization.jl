# Terminology

There are a lot of seemingly similar terms used for quantities in this package that refer to disparate things (or, are used slightly differently by other people in other places). Here, we try to best define these terms as we intend them.

## General Terms

### Feature

A quality or quantity associated with an atom that we wish to encode, such as atomic mass, row in the periodic table, etc.

### Encoding

The process of translating the value of a feature from its human-readable form (such as a float or a string) to whatever form will be ingested by a machine learning model. This could be as simple as an equality operation, but more often is, e.g. building a one-hot vector.

### Decoding

The inverse process to encoding. Note that in many cases (e.g. a continuous-valued feature encoded to a one-hot vector), the process isn't fully invertible, i.e. you can't get back a precise value but rather only a range corresponding to the associated onehot bin.

## Data types in `ChemistryFeaturization`

### Feature Descriptor

Describes the "features of a feature" – i.e. its name, possible values, instructions for encoding it, etc., but does NOT store an actual instance of its value. 

For more on the available types of feature descriptors, see [Feature Descriptors](@ref).

### Atoms Object

Describes a molecule, crystal, etc. in whatever representation will be ingested by an ML model (e.g. a graph), and can also store encoded features of that structure. 

For more on the available types of atoms objects, see [Atoms Objects](@ref).

### Featurization Object

Stores sets of feature descriptors and instructions for combining the values they encode on an atoms object.

For more on the available types of featurization objects, see [Featurizations](@ref).

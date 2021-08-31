# Primary Data Types

In [Terminology and Philosophy](@ref), we discussed general terms that pop up all over the package and the overall design philosophy of the package.

Let's go one step further and see some essential categories of types.

## Atoms

[Atoms](@ref) objects (terminology borrowed from [ASE](https://wiki.fysik.dtu.dk/ase/)) store information about the structure of a molecule, crystal, etc., and other useful information required to interpret the same.

## Feature Descriptor

[Feature Descriptors](@ref) represent the *"features of a feature"* - i.e., its name, possible values, instructions for encoding it, etc., but do NOT store an actual instance of their value.

They store all necessary information to encode and decode feature values on various parts of an [Atoms](@ref atoms) object and appropriately combine them into a single object (vector, matrix, etc.), describing the value/values of the feature for the entire object.

For example, if an [ElementFeatureDescriptor](@ref) encodes a vector for each atom in an object, they could be concatenated together into a matrix with a column for each atom to describe a structure.

Feature Descriptors must be designed with interoperability in mind. A `FeatureDescriptor` object must work deterministically with different datasets.

## Codec

[Codecs](@ref) define the actual encoding and decoding logic to be used with a feature descriptor.

## Featurization

A Featurization defines how the encoded values of different `FeatureDescriptors` can be collated and modified.

The implementation for a featurization must be as standardized and generic as possible, so that `FeatureDescriptors` that may be built using/use different datasets can be built and used easily.

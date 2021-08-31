# Primary Data Types

In [Terminology and Philosophy](@ref), we discussed general terms that pop up all over the package and the overall design philosophy of the package.

Let's go one step further and see some essential categories of types.

## Atoms

[Atoms](@ref) objects (terminology borrowed from [ASE](https://wiki.fysik.dtu.dk/ase/)) store information about the structure of a molecule, crystal, etc., and other useful information required to interpret the same.

## Feature Descriptor

[Feature Descriptors](@ref) represent the *"features of a feature"* - i.e., its name, possible values, instructions for encoding it, etc., but do NOT store an actual instance of their value.

## Codec

[Codecs](@ref) define the actual encoding and decoding logic to be used with a feature descriptor.

## Featurization

A [Featurization](@ref) stores sets of feature descriptors and defines the logic for combining the values they encode on an `Atoms` object.

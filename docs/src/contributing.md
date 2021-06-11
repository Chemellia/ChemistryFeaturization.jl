# Contributing

We welcome community contributions! This page contains some guidance to make the process go smoothly.

```@contents
Pages = ["contributing.md"]
```

## General Guidelines

For some high-level general guidance, please see [CONTRIBUTING.md](https://github.com/Chemellia/ChemistryFeaturization.jl/blob/main/CONTRIBUTING.md) on the GitHub repo.

The remainder of this page includes some pointers about package structure, etc. that could be helpful depending on what sorts of functionality you're interested in adding to the package. An understanding of "what goes where" will help in making sure you put your code in the right place! We also strongly encourage you to read the [Terminology/Philosophy](@ref) page, as it will help to understand the best ways to add different types of functionality.

TODO: flesh out everything below

## Implementing new atoms objects

* does it need to be a new type? Could it be accommodated by an existing type, or by an existing type with an expansion of functionality?
* make sure to subtype `AbstractAtoms`
* code belongs in `src/atoms/`, export statements belong in...
* make sure any/all sensible FD's and fzn's dispatch on it appropriately

## Implementing new feature descriptors

* does it need to be a new type? Could it be accommodated by an existing type, or by an existing type with an expansion of functionality?
* make sure to place it in the right place in the type hierarchy
* code belongs in `src/features/`
* make sure it encodes on as many atoms objects as are sensible

## Implementing new featurization schemes
* does it need to be a new type? Could it be accommodated by an existing type, or by an existing type with an expansion of functionality?
* should subtype `AbstractFeaturization`
* code belongs in `src/featurizations/`
* make sure it encodes on as many atoms objects as are sensible
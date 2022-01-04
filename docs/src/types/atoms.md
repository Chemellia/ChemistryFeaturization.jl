# [Atoms Objects](@id atoms)

```@docs
ChemistryFeaturization.Atoms
```

Atoms objects (terminology borrowed from [ASE](https://wiki.fysik.dtu.dk/ase/)) store information about the structure of a molecule, crystal, etc.

The parent abstract type is `AbstractAtoms`.

## AtomGraph

The `AtomGraph` type is used to store atomic graph representations. It can also be visualized using some customized formatting with the GraphPlot package.

```@docs
ChemistryFeaturization.Atoms.AtomGraph
ChemistryFeaturization.Atoms.AtomGraph(::String, ::String)
ChemistryFeaturization.Atoms.AtomGraph(::Crystal)
ChemistryFeaturization.Atoms.visualize
```

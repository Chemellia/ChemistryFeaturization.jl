# [Atoms Objects](@id atoms)

```@docs
Atoms
```

Atoms objects (terminology borrowed from [ASE](https://wiki.fysik.dtu.dk/ase/)) store information about the structure of a molecule, crystal, etc. as well as, optionally, encoded features and the featurization used to encode them. The parent abstract type is `AbstractAtoms`.

## AtomGraph
The `AtomGraph` type is used to store atomic graph representations. It can also be visualized using some customized formatting with the GraphPlot package.

```@docs
Atoms.AtomGraph
Atoms.AtomGraph(::String, ::String)
Atoms.visualize
```

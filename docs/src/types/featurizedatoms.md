# Featurized Atoms

A `FeaturizedAtoms` object is a flexible container for an `Atoms` object, a compatible `Featurization` object, and the resulting encoded features. The type is parameterized by the type of atoms object and the featurization type, e.g. `FeaturizedAtoms{AtomGraph,GraphNodeFeaturization}`.

This is intended to be the type that is directly fed into a Chemellia model, so the format of `encoded_features` should be whatever the associated model requires, and not necessarily human-readable (though interpretable through the `decode` function!).

```@docs
FeaturizedAtoms
```
# AtomFeat

The `AtomFeat` type stores featurization metadata and is carried along with feature vectors/matrices to ensure that encoded features (i.e. sequences of 0's and 1's) are always "decodable" to human-understandable quantities.

```@docs
ChemistryFeaturization.AtomFeat
```

## Constructors

```@docs
ChemistryFeaturization.AtomFeat(::Symbol, ::Bool, ::Integer, ::Real, ::Real, ::Bool)
ChemistryFeaturization.AtomFeat(::Symbol, ::Vector)
```

## Interfaces

```@docs
ChemistryFeaturization.build_featurization(::Vector{Symbol}; ::Vector{<:Integer}, ::Bool)
ChemistryFeaturization.make_feature_vectors(::Vector{AtomFeat})
ChemistryFeaturization.decode_feature_vector(::Vector, ::Vector{AtomFeat})
```

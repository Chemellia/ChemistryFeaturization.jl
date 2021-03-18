# AtomGraph

The `AtomGraph` type is used to store atomic graph representations.
It is a subtype of `LightGraphs.AbstractGraph{Float32}` (we've chosen to fix 32-bit floats for now to maintain compatibility with [Flux.jl](https://fluxml.ai/)), and has all the required functions defined such that most of the [LightGraphs](https://juliagraphs.org/LightGraphs.jl/latest/) API should work on these objects.

```@docs
ChemistryFeaturization.AtomGraph
```

## Constructors

```\
ChemistryFeaturization.AtomGraph(gr::SimpleWeightedGraph{Int32,Float32}, el_list::Vector{String})
ChemistryFeaturization.AtomGraph(gr::SimpleWeightedGraph{Int32,Float32}, el_list::Vector{String}, features::Matrix{Float32}, featurization::Vector{AtomFeat})
```

## Interfaces

```@docs
ChemistryFeaturization.add_features!
ChemistryFeaturization.add_features_batch!
ChemistryFeaturization.visualize_graph(ag::AtomGraph)
```

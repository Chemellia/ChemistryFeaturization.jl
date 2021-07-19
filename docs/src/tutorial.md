# Tutorial

```@contents
Pages = ["tutorial.md"]
Depth = 3
```

!!! note
    This tutorial is currently focused on [`AtomGraph`](@ref) objects, featurized using [`ElementFeatureDescriptor`](@ref)s in a [`GraphNodeFeaturization`](@ref), because these are the functionalities that are currently fully built, but we anticipate expanding it as more things are finished!

## Creating/reading in a structure
```@meta
DocTestSetup = quote
    using ChemistryFeaturization
end
```

### Create "from scratch"
We can build an `AtomGraph` "manually," by specifying an adjacency matrix and directly building the graph from that. Here I'll build a fictitious graph that's just an equiweighted triangle of carbon atoms:

```jldoctest
julia> adj_mat = Float32.([0 1 1; 1 0 1; 1 1 0]);

julia> triangle_C = AtomGraph(adj_mat, ["C", "C", "C"])
AtomGraph  with 3 nodes, 3 edges
   atoms: ["C", "C", "C"]

```
If you're working in an IDE that supports graphics output, you can also call `visualize(triangle_C)` to see the "ball-and-stick" graph.

### Reading from file
In a "real" application, you'll likely be reading structures from files such as .cif, .xyz, etc. Here, we'll read in the structure of WS<sub>2</sub>, downloaded from the [Materials Project](https://materialsproject.org):

```jldoctest WS2
WS2 = AtomGraph(joinpath("..", "files", "mp-224.cif"))

# output
AtomGraph mp-224 with 6 nodes, 9 edges
   atoms: ["W", "W", "S", "S", "S", "S"]

```
If you visualize this graph as above, you'll notice that it has two disconnected components. This isn't too surprising if we look at the 3D structure of this compound:
![WS2_structure](files/mp-224.png)
It's a two-dimensional material with two formula units per unit cell! Another way to see the disconnectedness of the graph is to index into the adjacency matrix in a particularly illustrative order:

```jldoctest WS2
WS2.graph[[1,4,6,2,3,5]].weights

# output
6×6 SparseArrays.SparseMatrixCSC{Float64, Int64} with 18 stored entries:
 1.0     0.9732   0.9732    ⋅       ⋅        ⋅ 
 0.9732  1.0      0.17143   ⋅       ⋅        ⋅ 
 0.9732  0.17143  1.0       ⋅       ⋅        ⋅ 
  ⋅       ⋅        ⋅       1.0     0.9732   0.9732
  ⋅       ⋅        ⋅       0.9732  1.0      0.17143
  ⋅       ⋅        ⋅       0.9732  0.17143  1.0
```

However, we have options in how we actually construct the graph.

## Building feature descriptors

## Building a featurization

## Featurizing structures

## Decoding encoded features
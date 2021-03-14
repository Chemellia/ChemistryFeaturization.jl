# Building Atomic Graphs

ChemistryFeaturization makes use of the `Atomic Simulation Environment` via Julia's `PyCall`. This allows it to build graphs from any file type that `ase.io.read()` can read in to an `Atoms` object – full list can be obtained [here](https://wiki.fysik.dtu.dk/ase/ase/io/io.html) or by running `ase info --formats` on the command line.

## Basic Options

The primary function that builds the graphs is (no surprise) `build_graph`.

```@docs
ChemistryFeaturization.build_graph
```

The primary one is the method of constructing the edge weights of the graph.\
The "cutoff" method corresponds to the method used in the [original cgcnn.py package](https://github.com/txie-93/cgcnn).
It considers all neighbors up to some cutoff radius (defaults to 8 Å), and adds up to a maximum of some number of neighbors (defaults to 12). In this implementation, the number cutoff is "soft," which is to say if there are additional neighbors at exactly the same distance, they will be added as well, even if the total number is above 12.\
The last option for this method is what function to use to set the edge weights based on separation distance. Currently the two options are `inverse_square` or `exp_decay`.

The other method is the Voronoi method (activated by setting `use_voronoi=true`), adopted from the [Ulissi Group's fork of cgcnn.py](https://github.com/ulissigroup/cgcnn).\
Note: that this method can only be used on fully periodic structures!\
It constructs a Voronoi tessellation to generate the neighbor list, and edge weights are added as the inverse of the polyhedra face areas.

## Batch Processing

The `build_graphs_batch` function allows for processing of sets of structure files with the same options outlined above.\

```@docs
ChemistryFeaturization.build_graphs_batch
```

Optionally, you can feed in a featurization scheme (vector of `AtomFeat` objects) and a set of feature vectors (dictionary from atomic symbols to prebuilt vectors, will be built automatically from the featurization scheme if not provided) if you'd like to featurize the graphs as well. It also takes all the same keyword arguments as the `build_graphs` function.

There is a corresponding function `read_graphs_batch(graph_folder::String)` which takes a path to a directory of serialized graphs (`graph_folder`) and reads them into an array of `AtomGraph` objects.

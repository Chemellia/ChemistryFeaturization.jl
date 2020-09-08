# ChemistryFeaturization.jl
Unified graph building and featurizing for Weave.jl, AtomicGraphNets.jl, and (maybe soon) more!

## Features

### Graph-building and featurization from CIF files
* Build graphs (as [SimpleWeightedGraphs](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl)) from CIF files using [PyCall](https://github.com/JuliaPy/PyCall.jl) to [pymatgen](https://pymatgen.org) functions
* Visualization using [GraphPlot](https://github.com/JuliaGraphs/GraphPlot.jl), check out the `visualize_graph` function in the `graph_functions.jl` file, you can make pretty pictures like these, whether the graph is simpler or more complicated (thickness of connections indicates weight of edge in graph (higher weights for nearer neighbors)):

<img src="img/graph_EuMgTl2.png" alt="graph_EuMgTl2" width="300" height="221"><img src="img/graph_K4W4O14.png" alt="graph_K4W414O14" width="305" height="221">

* Flexible featurization (currently onehot-style) and decoding: choose features to include, level of discretization, etc., and directly decode feature vectors to check values:
![](img/featurize_demo.gif)

### SMILES input
Sean to add...

## Requirements
* Julia 1.4+
* packages listed in `Project.toml`
* In addition, you will need your `PyCall` to have access to the `pymatgen` package, which can be added using `Conda.jl` as: `Conda.add("pymatgen"; channel="conda-forge")`, as well as the `rdkit` package (`Conda.add("rdkit"; channel="conda-forge")`)

## Future Plans:
### v0.1
* graph-building and featurizing for CrystalGraphConvNets.jl
* featurizing for Weave.jl
* featurization and inversion
* some basic documentation/examples

### later...
* saving out featurization metadata (features used, bins, etc.) in consistent format
* consistent handling of numerical vs. categorical features
* "hybrid" featurizations using features from multiple paradigms if available
* more input file formats? e.g. [SELFIES](https://github.com/aspuru-guzik-group/selfies)

# ChemistryFeaturization.jl
![Run tests](https://github.com/chemellia/ChemistryFeaturization.jl/workflows/Run%20tests/badge.svg)[![codecov](https://codecov.io/gh/chemellia/ChemistryFeaturization.jl/branch/main/graph/badge.svg?token=C0Fdt8BGnr)](https://codecov.io/gh/chemellia/ChemistryFeaturization.jl)

Flexible, modular, invertible featurization for Chemellia models including AtomicGraphnets.jl and WeaveModel.jl.

[Check out our JuliaCon 2021 talk here!](https://www.youtube.com/watch?v=la9asuZzVjU)

Documentation is starting to be built [here](https://chemellia.github.io/ChemistryFeaturization.jl/dev/), pardon the holes for the moment...

This package is currently focused on bulk systems. We have plans for more support of organic molecules, but at the moment, [MolecularGraph](https://github.com/mojaie/MolecularGraph.jl) is probably your best bet. PubChem stores many molecular features for the compounds they catalog, and their data can be accessed via [PubChemCrawler](https://github.com/JuliaHealth/PubChemCrawler.jl).

## Features

### Graph-building and featurization
* Build graphs from any file that Atomic Simulation Environment can read
* Visualization using [GraphPlot](https://github.com/JuliaGraphs/GraphPlot.jl), check out the `visualize_graph` function in the `graph_functions.jl` file, you can make pretty pictures like these, whether the graph is simpler or more complicated (thickness of connections indicates weight of edge in graph (higher weights for nearer neighbors)):

<img src="img/graph_EuMgTl2.png" alt="graph_EuMgTl2" width="300" height="221"><img src="img/graph_K4W4O14.png" alt="graph_K4W414O14" width="305" height="221">

* Flexible featurization (currently onehot-style) and decoding: choose features to include, level of discretization, etc., and directly decode feature vectors to check values

* Much more coming!

## Requirements
* Julia >= 1.4
* packages listed in `Project.toml`
* In addition, you will need your `PyCall` to have access to the `pymatgen` package, which can be added using `Conda.jl` as: `Conda.add("pymatgen"; channel="conda-forge")`, as well as the `rdkit` package (`Conda.add("rdkit"; channel="conda-forge")`)



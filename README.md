# ChemistryFeaturization.jl
![Run tests](https://github.com/chemellia/ChemistryFeaturization.jl/workflows/Run%20tests/badge.svg) [![codecov](https://codecov.io/gh/chemellia/ChemistryFeaturization.jl/branch/main/graph/badge.svg?token=C0Fdt8BGnr)](https://codecov.io/gh/chemellia/ChemistryFeaturization.jl)  [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://chemistryfeaturization.chemellia.org/stable/)
 [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://chemistryfeaturization.chemellia.org/dev/)

**This branch has an example used in some seminar demonstrations and was made so that the exact version of the code used would be easily accessible. If you are here for that code, the relevant file is `examples/seminar_demo.jl`.**

Interface for modular, flexible, invertible featurization of atomic structures for machine learning purposes.

[Check out our JuliaCon 2021 talk here!](https://www.youtube.com/watch?v=la9asuZzVjU) (note that some syntax is out of date)

This package is currently focused on bulk systems. We have plans for more extensive support of organic molecules, but at the moment, [MolecularGraph](https://github.com/mojaie/MolecularGraph.jl) is probably your best bet. PubChem stores many molecular features for the compounds they catalog, and their data can be accessed via [PubChemCrawler](https://github.com/JuliaHealth/PubChemCrawler.jl).

## Features

For a fairly extensive tour, check out the tutorial in the docs!


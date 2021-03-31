# ChemistryFeaturization.jl

*Data types and featurization schemes for Alchemy models and beyoond!*

ChemistryFeaturization.jl is meant to be a unified interface for translating atomic structures (molecules, crystals, etc.) into data structures and sets of features to be used in machine learning models provided by the Alchemy suite of packages, such as [AtomicGraphNets.jl](https://github.com/aced-differentiate/AtomicGraphNets.jl).

This package is in development as part of the [ACED project](https://www.cmu.edu/aced/), funded by ARPA-E DIFFERENTIATE and coordinated by [Carnegie Mellon University](https://www.cmu.edu/), in collaboration with [Julia Computing](https://juliacomputing.com/), [Citrine Informatics](https://citrine.io/), and [MIT](https://web.mit.edu/).

## Purpose

This is intended to serve as a "helper package" of sorts for other packages that actually build the models, and NOT as a standalone tool (for instance, ChemistryFeaturization.jl it does not itself implement any models).

It provides...

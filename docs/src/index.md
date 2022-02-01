# ChemistryFeaturization.jl

ChemistryFeaturization.jl provides an interface for translating atomic structures (molecules, crystals, etc.) into data structures and sets of features to be used in machine learning models, in particular those provided by the Chemellia suite of packages, such as [AtomicGraphNets.jl](https://github.com/Chemellia/AtomicGraphNets.jl).

To learn more about the package, we would suggest starting by reading this page for a high-level perspective, and then proceeding to the tutorials. Of course, full docstrings are also available at [API](@ref), likely useful particularly if you want to implement new feature descriptor and/or featurization types for your own models.

## Design Philosophy
Several aspects of the design of ChemistryFeaturization make it stand out:
1. **Flexibility**: Other packages for featurizing atomic structures can be rigid, demanding a specific set of features be encoded in a particular way. In ChemistryFeaturization, there may be some convenience function dispatches for defaults, but you always have the ability to tweak which features you include, how they're encoded, etc.
2. **Decodability**: Typical encoding schemes are opaque regarding what exactly has been encoded and with what fidelity. The structure of ChemistryFeaturization is such that there is always the ability to decode any encoded value to inspect.
3. **Extensibility**: Fundamentally, ChemistryFeaturization is an interface. It includes a small number of concrete implementations that we envision will be of broad utility, but its power will really come from how it facilitates sharing and reuse while maintaining the two characteristics summarized above.

In addition, as the [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) interface grows and hopefully gains a user base, more integration with that ecosystem is planned.

## API Overview
This section is meant to give a high-level summary of the types and functions defined in ChemistryFeaturization. For details on syntax and configuration options, please see the associated pages/docstrings.

### Feature Descriptors
A **feature descriptor** is a core unit describing a single feature. It could be a feature of a single atom (e.g. atomic mass, electronegativity, ...), a pair of atoms (e.g. bond length), larger collections such as a local environment, or an entire structure. All feature descriptors are subtypes of [`AbstractFeatureDescriptor`](@ref), and need to dispatch the [`get_value(feature, structure)`](@ref) function to compute the value(s) of the feature for a given structure, as well as ideally a few others such as [`encodable_elements`](@ref) and [`default_codec`](@ref). An example feature descriptor type (and, in fact, the only one with a concrete implementation in ChemistryFeaturization) is `ElementFeatureDescriptor`, a subtype of `AbstractAtomFeatureDescriptor` that uses a lookup table to compute feature values based only on atomic symbol. Currently, other feature descriptors are implemented in other packages that use the ChemistryFeaturization interface.

### Codecs
A **codec** (encoder/decoder) describes how values of a feature should be encoded, and also allows for them to be decoded. Typically, one feature descriptor is paired with one codec in a featurization scheme. All codecs are subtypes of `AbstractCodec` and need to dispatch [`encode(val, codec)`](@ref) and [`decode(val, codec)`](@ref). Examples of codecs include `OneHotOneCold` (for standard bitstring encodings of categorical features, and with flexible binning schemes for continuous-valued ones) and `DirectCodec` (only compatible with numerical-valued features and simply scales their value by some constant).

### Featurizations
A **featurization** (or featurization scheme) is a collection of feature descriptors and associated codecs, as well as a way to combine the outputs of those codecs into whatever shape/format is needed for ingestion by a given model. Typically, we anticipate a one-to-one pairing between featurizations and model architectures. For example, [AtomicGraphNets.jl](https://github.com/Chemellia/AtomicGraphNets.jl) uses the `GraphNodeFeaturization`, which encodes and concatenates a set of atom features for every node in a graph, and "stacks" the resulting vectors into a matrix, which, along with the adjacency matrix of the graph, serves as the input into a crystal graph model.

All featurizations are subtypes of `AbstractFeaturization` and should dispatch the [`features`](@ref) function to return a list of feature descriptors. They may also need to dispatch [`encode`](@ref) and [`decode`](@ref) if the default dispatches are not the desired behavior.

### `FeaturizedAtoms` objects
A `FeaturizedAtoms` object is just a container for a representation of an atomic structure, a featurization, and the resulting encoded features from that featurization applied to the atomic structure. It is the return type of `featurize(atoms, featurization)` and has two type parameters for the type of the structure representation and the type of featurization.

`FeaturizedAtoms` objects are useful for prefeaturizing and storing/serializing data for feeding into a model later, while ensuring all the necessary metadata for decoding is still attached.

## Acknowledgements
This package is in development as part of the [ACED project](https://www.cmu.edu/aced/), funded by ARPA-E DIFFERENTIATE and coordinated by a team from [Carnegie Mellon University](https://www.cmu.edu/), in collaboration with [Julia Computing](https://juliacomputing.com/), [Citrine Informatics](https://citrine.io/), and [MIT](https://web.mit.edu/). [Dr. Rachel Kurchin](https://rkurchin.github.io) is the lead developer, and is grateful to [MolSSI](https://molssi.org) for support of the first six months of development. [Anant Thazhemadam](https://thazhemadam.github.io/blog/) and [Dhairya Gandhi](https://github.com/DhairyaLGandhi) have also been major contributors of both code and ideas.

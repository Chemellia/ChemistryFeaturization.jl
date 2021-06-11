# Changelog

I'm generally trying to adhere to [semver](https://semver.org) here. This means that in v0.*, assume breaking changes are always possible, even without a major version bump...however, I will try to always note them here if they happen...

Categories to include for each release, if relevant: breaking, added, fixed, removed/deprecated

## v0.3.0
Major restructure to make extensibility far easier and allow better sharing of elements between featurization schemes! Also streamlined the codebase a lot by organizing things into sensible modules, so technically basically everything is a breaking change because all of the imports are from the submodules now.

### Breaking/Added
* create modules for different functionalities: `Atoms`, `Features`, `Featurizations`, `Utils`, etc.
* separate out encoding/decoding functionality into a new `Codec` type and associated `Codecs` module
* merged functionality of `AtomGraph` and `WeaveMol` into `AtomGraph` object, with a more generic `encoded_features` field

## v0.2.2 [2021-02-22]

### Added
* updated pretty printing for `AtomGraph` to include `id` field

### Fixed
* added some docstrings, updated others

## v0.2.1 [2021-02-22]
### Breaking
* rename `build_atom_feats` to `build_featurization` to avoid ambiguity with `make_feature_vectors`

### Fixed
* proper docstrings for more things incl. `AtomFeat`, `AtomGraph`, `add_features!`

## v0.2.0 [2021-02-16]
### Breaking
* add `id` field to `AtomGraph` (should only be breaking for reading in serialized graphs)
* `build_graphs_batch` now returns list of graphs and optionally serializes rather than always serializing and never returning

### Added
* `add_features_batch!` function
* `AtomGraph` constructor directly from adjacency matrix (previously was only from a `SimpleWeightedGraph`)

### Fixed
* remove deprecated syntax of opening DataFrames from CSV
* made separate files into modules to avoid redefinition warnings during precompilation

## v0.1.1 [2021-02-10]

### Fixed
* add check for NaN values in graph laplacian

## v0.1.0 [2020-12-22]
Initial release!

### Added
* Create AtomGraph and AtomFeat types
* Basic graph visualization functions
* Graph-building from CIF files via the cgcnn.py "cutoff" method, with support for nonperiodic systems as well

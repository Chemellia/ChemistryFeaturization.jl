# Changelog

I'm generally trying to adhere to [semver](https://semver.org) here. This means that in v0.*, assume breaking changes are always possible, even without a major version bump...however, I will try to always note them here if they happen...

Categories to include for each release, if relevant: breaking, added, fixed, removed/deprecated

## Upcoming
### Added
* add pretty printing for `GraphNodeFeaturization`, `FeaturizedAtoms`, uniformize spacing via `\t` across these and `AtomGraph`
* [logic] move `data/` directory to root of repo
* [logic] move `build_onehot_vec` to `OneHotOneCold`

### Fixed
* export `featurize` and `decode` properly for `FeaturizedAtoms`

## v0.4.2 [2021-07-02]
### Added
* export `featurize` and `decode` functions at top-level module

### Fixed
* `id` is positional rather than keyword argument in `AtomGraph` constructor so that broadcast works properly, also it defaults to filename sans extension when constructing from file

## v0.4.1 [2021-07-01]
### Fixed
* add back wayward semicolon so that keyword arguments work in `AtomGraph` constructor from file

## v0.4.0 [2021-06-30]
### Added
* create `FeaturizedAtoms` type, remove `featurization`, `encoded_features` fields from Atoms objects
* add docstrings for various things, add section of docs for Codecs

### Fixed
* macOS CI fixed, but possibly not in an ideal way because it seems to result in occasional precompile warnings upon update
* Docs for stable/tagged versions now build properly

### Removed/Deprecated
* remove encoded features field from `AtomGraph` type

## v0.3.1 [2021-06-18]
### Fixed
* import from submodules before top-level export so that exports actually work, d'oh
* remove broken dimension check in `AtomGraph` constructor (doesn't work with generic featurization)
* fix type assertion of `nbins` in `GraphNodeFeaturization` constructor (`Vector{<:Integer}`, NOT `Vector{Integer}`)

## v0.3.0 [2021-06-17]
Major restructure to make extensibility far easier and allow better sharing of elements between featurization schemes! Also streamlined the codebase a lot by organizing things into sensible modules, so technically basically everything is a breaking change because all of the imports are from the submodules now.

### Breaking/Added
* create modules for different functionalities: `Atoms`, `Features`, `Featurization`, `Utils`, etc.
* separate out encoding/decoding functionality into a new `AbstractCodec` type and associated `Codec` module
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

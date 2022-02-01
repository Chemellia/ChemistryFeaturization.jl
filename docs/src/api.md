# API
```@contents
Pages = ["api.md"]
Depth = 5
```

## Multiple Dispatch
By design, several functions have dispatches on multiple abstract types. Their docstrings are presented here for an overview, and then specific cases are duplicated in relevant sections below.

```@docs
encode
decode
encodable_elements
output_shape
```

## Feature Descriptors
### Abstract Types
```@docs
AbstractFeatureDescriptor
AbstractAtomFeatureDescriptor
AbstractPairFeatureDescriptor
```

### Working with Feature Descriptors
```@docs
get_value
encodable_elements(::AbstractFeatureDescriptor)
default_codec
encode(atoms, ::AbstractFeatureDescriptor)
decode(atoms, ::AbstractFeatureDescriptor)
```

### `ElementFeature` submodule
The `ElementFeature` submodule includes a concrete implementation of an `AbstractAtomFeatureDescriptor` for the case of features whose values can be computed from a lookup table of elemental symbols. It also has some utility functions for defining default encoder settings.

```@docs
ElementFeature.ElementFeatureDescriptor
ElementFeature.fea_minmax
ElementFeature.default_log
ElementFeature.default_categorical
ElementFeature.get_bins
ElementFeature.get_param_vec
```

To see the concrete implementations of interface functions such as `get_value`, `encodable_elements`, and `default_codec` on `ElementFeatureDescriptor` objects, take a look at the source code at `src/features/elementfeature.jl`.

The `ElementFeature` module also makes use of the `Data` module, which provides data to populate lookup tables for a variety of commonly-desired features. In particular, it provides the constants `element_data_df` which will be automatically used as the lookup table for an `ElementFeatureDescriptor` if none is provided, as well as `elementfeature_info`, a dictionary providing information about the available features included in `element_data_df`. 

## Codecs
### Codec Interface
```@docs
AbstractCodec
encode(val, ::AbstractCodec)
decode(val, ::AbstractCodec)
output_shape(::AbstractCodec)
```

### Concrete Types
The `OneHotOneCold` codec is a very common encoding scheme. For categorical-valued features, it encodes a bitstring with a length equal to the number of possible values composed of zeros except in the slot corresponding to the value in that instance. For continuous-valued features, a binning scheme must be specified, or defaults will be chosen using the utility functions shown below.
```@docs
OneHotOneCold
DirectCodec
```

## Featurizations

```@docs
AbstractFeaturization
features
encodable_elements(::AbstractFeaturization)
encode(atoms, ::AbstractFeaturization)
decode(atoms, ::AbstractFeaturization)
```

## `FeaturizedAtoms` objects
```@docs
FeaturizedAtoms
featurize
decode(::FeaturizedAtoms)
```

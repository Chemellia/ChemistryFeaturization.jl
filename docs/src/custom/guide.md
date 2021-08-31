# New types

A quick and nifty guide to getting started with writing new types.

* It is recommended to have custom pretty printing for types - especially Atoms, FeatureDescriptors, and Featurizations

## Atoms Objects

* Always be a descendant `AbstractAtoms`

## Codecs

* Always be a descendant `AbstractCodec`
* two fields - `encode_f::Function` and `decode_f::Function`

Remember, a Codec is defined by the parameters the encoding/decoding scheme it is built to support takes.

## Feature Descriptors

* Always be a descendant `AbstractFeatureDescriptor`
* have a field of type Codec (`::AbstractCodec`)
* an appropriate `encodable_elements`
* callable syntax that can be dispatched generically over any AbstractAtoms object
* `encode`, and `decode`
* `output_shape`
* A default Codec
* Callable syntax for the default codec, dispatched with the FeatureDescriptor type as an argument
* default encode-decode functions for compatible Codec(s)

## Featurizations

* Always be a descendant `AbstractFeaturization`
* `encode` and `decode`
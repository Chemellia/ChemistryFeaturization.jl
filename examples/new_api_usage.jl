# for features with built-in data, there will be convenience constructors that can build AtomFeats with default binning, etc. just from the names, so you could do this...
feature_names = ["Group", "Row", "Block", "Atomic mass", "Atomic radius", "X"]
features = AtomFeat.(feature_names)

# or, of course, have a higher degree of control by specifying bin numbers, etc. and then just doing something like...
features = AtomFeat.(zip(feature_names, num_bins, logspaced))

# the point is, rarely should the user ever have to write the encode_f and decode_f themselves (really only if they're defining their own custom features)

# then you could build a featurization like this (default `combine` would just be some concatenation)
featurization = GraphNodeFeaturization(features)

# or even build it directly from the feature names and also be able to pass in the custom binnings, etc. (this becomes the analogue of `build_featurization` from before)
# and we no longer need `make_feature_vectors` because the GraphNodeFeaturization builds and stores them when it's constructed
featurization = GraphNodeFeaturization(feature_names)

# then read in the graphs (assume we're not constructing them here because the mechanics of that won't change too much anyway) and featurize them...
graph_paths = do_something.(info_df)
inputs = deserialize.(graph_paths)
featurize!(inputs, featurization)
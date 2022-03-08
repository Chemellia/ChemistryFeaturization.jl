## Featurizing some AtomGraphs (this is mostly a compressed version of https://chemistryfeaturization.chemellia.org/dev/tutorial/atomgraphs/, so comments here are minimal as they would just replicate those)
# (don't forget to activate the environment for this repo! I've committed the manifest on this branch so you should be able to duplicate mine exactly)

using AtomGraphs
using AtomicGraphNets
using ChemistryFeaturization
using ChemistryFeaturization.ElementFeature

adj_mat = Float32.([0 1 1; 1 0 1; 1 1 0])
triangle_C = AtomGraph(adj_mat, ["C", "C", "C"])
visualize(triangle_C)

WS2 = AtomGraph(pathof(ChemistryFeaturization)[1:end-29]*"test/test_data/strucs/mp-224.cif")
visualize(WS2)

# categorical feature
block = ElementFeatureDescriptor("Block")
block(triangle_C) # equivalent to get_value(block, WS2)
block_codec = default_codec(block)
d_encoded = encode("d", block_codec)
decode(d_encded, block_codec)
encode("s", block_codec)
WS2_block_encoded = encode(WS2, block)
decode(WS2_block_encoded, block)

# continuous-valued feature
amass = ElementFeatureDescriptor("Atomic mass")
amass(triangle_C)
triangle_C_amass_encoded = encode(triangle_C, amass)
decode(triangle_C_amass_encoded, amass)
amass_codec_hires = OneHotOneCold(false, get_bins(amass_codec.bins, nbins=16, logspaced=true))
decode(encode(triangle_C, amass, amass_codec_hires), amass_codec_hires)

# custom feature
using DataFrames
lookup_table = DataFrame(["C" 42; "As" 0], [:Symbol, :MeaningOfLife]); # make a custom lookup table for another feature
meaning = ElementFeatureDescriptor("MeaningOfLife", lookup_table)
meaning(triangle_C)
meaning(WS2)

# building a featurization
fzn = GraphNodeFeaturization(["Block", "Atomic mass", "X"])
featurized_WS2 = featurize(WS2, fzn)
propertynames(featurized_WS2)
decode(featurized_WS2)

## Implementing the interface for an AtomsBase structure
using AtomsBase
using AtomIO

struct AtomsBaseFeaturization
    features::Vector{<:AbstractAtomFeatureDescriptor}
    codecs::Vector{<:AbstractCodec}
end

AtomsBaseFeaturization(gnf::GraphNodeFeaturization) = AtomsBaseFeaturization(gnf.features, gnf.codecs)

features(abf::AtomsBaseFeaturization) = abf.features

atoms = load_system(pathof(ChemistryFeaturization)[1:end-29]*"test/test_data/strucs/mp-224.cif")

encode(atoms, fzn)
### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 8cc13d04-83a9-11ec-2d38-47d50cf23b83
md"# ChemistryFeaturization Usage Example: Featurizing `AtomGraph`s
This tutorial will demonstrate usage of ChemistryFeaturization.jl for a case that has already implemented the interface. The structure representation is `AtomGraph`, provided by the [AtomGraphs](https://github.com/Chemellia/AtomGraphs.jl) package, and the featurization is `GraphNodeFeaturization` from [AtomicGraphNets](https://github.com/Chemellia/AtomicGraphNets.jl). We will choose some `ElementFeatureDescriptors` by which to featurize each node of our `AtomGraph` object and step through how to do it.
"

# ╔═╡ 50fbffca-95d9-4b2b-955e-b5453574ddaf
using ChemistryFeaturization, AtomGraphs, AtomicGraphNets, PlutoUI

# ╔═╡ db3022ee-9241-4446-87ad-92c8359d1e31
#! format: off
TableOfContents()

# ╔═╡ 96095ce5-6edb-46f0-9884-397b318ed631
md"## Creating `AtomGraph` objects
`AtomGraph` objects are simply adjacency matrices with elemental symbols labeling each node. We can create one \"from scratch\" by manually specifying an adjancey matrix, like so...
"

# ╔═╡ 0f71df51-9e39-4768-8965-d79c0105f3e8
adj_mat = Float32.([0 1 1; 1 0 1; 1 1 0])

# ╔═╡ 84d37c16-1afc-4ad9-9ae5-8beb5c8d3bc8
triangle_C = AtomGraph(adj_mat, ["C", "C", "C"])

# ╔═╡ b2263a81-f2f0-4b0a-a6f0-393e36c1e0ba
md"The built-in `visualize` function in AtomGraphs.jl doesn't work in Pluto notebooks, so we'll create a slightly modified one that will display inline properly and visualize the graph.
"

# ╔═╡ 8add5136-1a19-440c-a2d7-0d19f2d5e809
begin
	using GraphPlot, Graphs
	
	function plutoviz(ag::AtomGraph)
		sg = SimpleGraph(ag.graph.weights)
		gplot(
        	sg,
        	nodefillc = AtomGraphs.graph_colors(elements(ag)),
        	nodelabel = elements(ag),
        	edgelinewidth = AtomGraphs.graph_edgewidths(ag),
    	)
	end
end

# ╔═╡ fc212e05-bc3b-484c-84e4-e051f7c555df
plutoviz(triangle_C)

# ╔═╡ 3fb29959-a7db-4e48-a527-e81ac85b486d
md"Of course, in practice we'll more likely be reading in structures and building graphs from files such as .cif, .xyz, etc. Here, we'll read in the structure of WS<sub>2</sub>, downloaded from the [Materials Project](https://materialsproject.org):
"

# ╔═╡ 718c6c24-05ff-4aca-bb54-787992e1a517
WS2 = AtomGraph("../files/mp-224.cif")

# ╔═╡ a501a4db-fed5-43e7-bbc2-bf5657475ee1
md"The graph is automatically assigned an `id` based on the filename it was read from, but you can pass a value to override this and name it something else.

We can visualize it just like we did before:"

# ╔═╡ 114a2eff-acef-42b1-bedc-96dfc260e127
WS2

# ╔═╡ 084f5ad9-6381-4663-97f6-f5f73ac2944d
plutoviz(WS2)

# ╔═╡ ee8549d0-6ba2-4749-be52-74814756deca
md"Interesting! It seems this graph has two disconnected components. This isn't too surprising if we look at the 3D structure of this compound:"

# ╔═╡ 6af1e757-79ce-4e5c-88f1-45b4d39a993f
begin
	using ImageShow, Images
	load("../files/mp-224.png")
end

# ╔═╡ 75b26589-b74f-4508-8d62-fbf40423d64a
md"It's a two-dimensional material with two formula units per unit cell! Another way to see the disconnectedness of the graph is to index into the adjacency matrix in a particularly illustrative order:"

# ╔═╡ 39e9b99a-3964-4f5d-a5a7-4bf3726c4429
WS2.graph[[1,4,6,2,3,5]].weights

# ╔═╡ 9ffc63ac-04cd-4524-88ad-2339035c8daf
md"However, we have options in how we actually construct the graph. The default option is based on the scheme from [the original cgcnn.py implementation](https://github.com/txie-93/cgcnn), which essentially involves setting a maximum neighbor distance and a maximum number of neighbors. However, in contrast to that implementation, we construct weighted graphs (with the user having an ability to specify the weight decay function with separation distance; it defaults to inverse-square).

An arguably more physical way to construct neighbor lists and graphs is by doing a [Voronoi partition](https://en.wikipedia.org/wiki/Voronoi_diagram) of the atomic coordinates. In this scheme, the neighbor list of an atom is any atom with which its Voronoi polyhedron shares a face, and the edge weights can be determined using the areas of the faces. Let's try that with our WS<sub>2</sub> structure..."

# ╔═╡ 12e7461c-af53-4449-a800-4663fb4839ee
WS2_v = AtomGraph(joinpath("..", "files", "mp-224.cif"), use_voronoi=true)

# ╔═╡ 75100b89-f2cc-40bc-9d8a-36f52831292e
WS2_v.graph[[1,4,6,2,3,5]].weights

# ╔═╡ a9be0242-65ff-465d-8c1e-6e3792e38147
md"### Batch Processing
One final note for this section: 
the `AtomGraph` constructor broadcasts! So if you have a directory full of structure files (say, `strucs/`), you can get a list of `AtomGraph` objects by:
```julia
ags = AtomGraph.(readdir(\"strucs/\", join=true))
```
"

# ╔═╡ 0ee1995b-77a8-4415-b686-f876fa2761f8
md"## Building and Encoding Feature Descriptors
What types of features of our structure do we want to encode in our graph? Let's keep things simple for now and consider features that can be encoded only by knowing the elemental identity of a given atom (node in our graph). The package includes a bunch of built-in data, and you can also provide your own for features we haven't included!

We can easily construct these for built-in features...
"

# ╔═╡ ee15cd74-bdc2-4437-a100-1d5857f33b7e
using ChemistryFeaturization.ElementFeature

# ╔═╡ 4f73e95c-ffdd-4227-bbe5-14b80413cd46
md" ### Categorical features
Let's start with a categorical feature, that is, one that takes on a finite set of discrete values. One example of this is which block (s, p, d, or f) in the periodic table an element resides in."

# ╔═╡ 773c4a97-c5ab-486e-8286-08b3b7167e7d
block = ElementFeatureDescriptor("Block") # categorical feature denoting s-, p-, d-, or f-block elements

# ╔═╡ 26eb23af-d1b2-43be-8f04-f820f441c8e3
md"We can get the values of the feature for a given structure by \"calling\" it directly..."

# ╔═╡ 25909b37-9a83-44a3-a287-a19c7e542391
block(triangle_C)

# ╔═╡ 22fe0263-ffc3-4cbe-9f7b-29f3a0ce6d97
md"...or by using the `get_value` function."

# ╔═╡ da01488b-efd0-4653-a9a2-b8defe93ff17
get_value(block, WS2)

# ╔═╡ bb6bcb31-2905-49af-8070-bb493069cc84
md" Of course, vectors of single characters are not going to be all that useful to feed into a machine learning model. To \"translate\" these human-readable values, we need to _encode_ them. For this, we'll use a codec object. Codec is short for \"encoder-decoder\", reflecting a key design principle of ChemistryFeaturization that we should always know what information we've encoded and be able to invert that encoding process.

A common method of encoding categorical-valued featuers is using so-called \"one-hot\" encoding. In ChemistryFeaturization, this is implemented via the `OneHotOneCold` codec. We can retrieve a \"sensible\" one for our `block` feature using the `default_codec` function...
"

# ╔═╡ 644b1a5a-65e7-4731-8fd8-058670636b1d
block_codec = default_codec(block)

# ╔═╡ 2ac9fd2a-44c6-4c40-b8dc-9ce47c8eb333
md"So what is this `OneHotOneCold` thing? It may be easiest to see by using it. We can encode single values..."

# ╔═╡ 1127b60a-b52d-4d30-8d4a-06a38dba4d04
d_encoded = encode("d", block_codec)

# ╔═╡ d6205d96-6cb3-438b-8adb-b5d6b8acaa45
s_encoded = encode("s", block_codec)

# ╔═╡ 32fa1fd9-55cd-4615-a2bc-2a9e0b2524c1
md"...or even whole structures..."

# ╔═╡ 47cd8c6a-32ba-4f23-96bf-8b678c8b1721
WS2_block_encoded = encode(WS2, block)

# ╔═╡ 99a39408-d1f8-4afb-8f18-7bcd1b7faf08
md"The output for a single value is a bitstring of 0's with a 1 in the \"slot\" corresponding to the value, where the \"slots\" are specified by the `vals` field of the codec, in this case `[\"d\", \"f\", \"p\", \"s\"]`. For a structure like our `WS2 AtomGraph`, these vectors are concatenated into a matrix where the first index is the index of the atom and the second indexes into the feature vector. Note that we could have called `encode(WS2, block, block_codec)` above for the same result, but `encode` will call `default_codec` by default so it's not necessary.

Note that we can always `decode` what we `encode`...(and ChemistryFeaturization will also internally call `default_codec` so if you're using the defaults, you can feed either the codec or the feature descriptor)"

# ╔═╡ f9d92884-cdce-4b9c-984f-3f5bd69ea020
decode(s_encoded, block_codec)

# ╔═╡ 8ed6f27e-d40f-4d05-a92c-3515980aad59
decode(WS2_block_encoded, block)

# ╔═╡ bbcca9da-3027-4025-885f-21b12454a5aa
md"So what's that other field in the codec that has a value of `true`? Read on..."

# ╔═╡ b00a2447-06bc-4d93-b717-96f01771c176
md" ### Continuous-valued features
Some features, such as atomic mass, are better described as continua of values. One example might be the mass of an atom, another built-in `ElementFeatureDescriptor`:
"

# ╔═╡ b94b89c1-5496-41d5-be02-6136763f1153
amass = ElementFeatureDescriptor("Atomic mass") # continuous-valued feature

# ╔═╡ b1595f7e-21c1-44d7-b942-2d652cb5a28b
amass_codec = default_codec(amass)

# ╔═╡ 51db299f-5a1c-431c-a248-42be75a1e18f
md"Note that this time that first flag is `false`. The flag describes whether the codec is describing a categorical or continuous-valued feature, and this influences how the `bins` are interpreted. For categorical features, each value in `bins` corresponds to a possible value of the feature. For continuous features, `bins` represents N+1 *edges* of N bins into which we've divided the possible values of this feature. 

In addition, in this case of atomic mass, the bins are logarithmically spaced by default. For more on how to tune these defaults, check out the `onehotonecold_utils.jl` source file.

Note a further consequence of this difference in the interpretation of `bins` for categorical vs. continuous features:"

# ╔═╡ e436c12d-c5fe-4904-be9c-6802559a5fd1
length.([block_codec.bins, amass_codec.bins])

# ╔═╡ 1093a1fb-edfb-4176-b310-4e1fa5572473
output_shape.([block_codec, amass_codec])

# ╔═╡ 6d75ebb9-3212-4058-8236-8c900b961019
md"(Beware OBO errors all ye who enter here.)

Let's try encoding/decoding some values!"

# ╔═╡ 92c9954b-72e7-4d91-988a-915b2c5999d5
amass(triangle_C)

# ╔═╡ 91afd9c2-fc59-4b2b-b1f0-b9cf61a06ade
triangle_C_amass_encoded = encode(triangle_C, amass)

# ╔═╡ 2c849e35-9161-4ad4-aa80-6ee203205810
decode(triangle_C_amass_encoded, amass)

# ╔═╡ c4a25ad0-ee4f-43f9-9040-5d1d39aaa72a
md"Oooh, interesting, what's happened here?

Turns out one-hot encoding loses some information for continuous-valued features, so instead of getting back carbon's atomic mass of 12.01, we can back tuples indicating the edges of the bins. If you want to encode with higher resolution, you can build a codec with more bins!"

# ╔═╡ 317cd6d0-18b7-4fcd-8137-46d6ed2bdbff
amass_codec_hires = OneHotOneCold(false, get_bins(amass_codec.bins, nbins=16, logspaced=true))

# ╔═╡ c307dc74-879e-45d7-8cd1-076580a66c43
decode(encode(triangle_C, amass, amass_codec_hires), amass_codec_hires)

# ╔═╡ 6e61b457-12dd-4038-8e4c-23bc3cc62ba5
md" ### Custom features"

# ╔═╡ 8ccfc9bc-0b4b-4ea7-89c7-86b7d3b756e4
md"But suppose you have another feature that's not included. You can easily provide a lookup table (or even an entire custom encoding function!) yourself, like so..."

# ╔═╡ c10632e4-4f50-42a1-93e7-e4736ab0237e
begin
	using DataFrames
	lookup_table = DataFrame(["C" 42; "As" 0], [:Symbol, :MeaningOfLife]); # make a custom lookup table for another feature
	meaning = ElementFeatureDescriptor("MeaningOfLife", lookup_table)
end

# ╔═╡ 51d2f6d9-681e-46a4-89d8-d0d25a3fdebd
md"Once you've done this, you can use this feature the same way you would use the built-in ones. Just make sure you're okay with the default decisions the package makes about whether to treat the feature as categorical vs. continuous and space the bins linearly or logarithmically (again, see `src/codecs/onehotonecold_utils.jl` for details) and tweak the flags if you're not."

# ╔═╡ 165f0c92-0f01-4525-a678-51178f36c328
md"## Building a featurization
Next, we can combine these feature descriptors into a _featurization object_, which allows convenient encoding of multiple features on a structure, and also combining of those encoded features in a manner appropriate for feeding into a model. In the case of `GraphNodeFeaturization`, we construct a vector for each node in an `AtomGraph` by concatenating encoded features together, and then stack these vectors to form a feature matrix that we could feed into an AtomicGraphNets model.

This featurization has a convenience constructor that will build the `ElementFeatureDescriptor`s if you just pass in names of features, but with our custom lookup table feature, we would need to construct it by directly passing the feature descriptors:
"

# ╔═╡ 498144b5-e599-45da-97b7-a65931429ab6
fzn = GraphNodeFeaturization([block, amass, meaning])

# ╔═╡ e8d10af1-bde0-49ee-ac9a-33c82da7cd1c
md"(As a quick side note, the featurization is basically just a bundle for feature descriptors and associated codecs, which we can of course inspect:)"

# ╔═╡ 9faaa0b2-ef70-470f-9011-83a1367e6963
fzn.codecs

# ╔═╡ 44a3a9bb-7c6b-4270-bea8-71ecb3fa7aed
md"## Featurizing structures
We've already seen how we can retrieve and encode values of individual features for atomic structures. We can of course do this for featurizations too:
"

# ╔═╡ 96bb6f6d-ada0-4fb8-bdf2-f12de1e68b71
encode(triangle_C, fzn)

# ╔═╡ 2a619e89-8398-48e7-afb3-ba11a9557820
md"
If we want to attach the encoded features to the graph, we can use the `featurize` function, which returns a `FeaturizedAtoms` object. This is the recommended approach for preprocessing many structures, as a serialized `FeaturizedAtoms` object stores the structure, the featurization, and the encoded features so you know what you've encoded it and how, and have the ability to decode.
"

# ╔═╡ 06c294c4-fa86-4d17-891d-39af1c81586e
try # hide
	featurize(WS2, fzn)
catch e # hide
	@show e # hide
end # hide

# ╔═╡ d5535021-8c23-4640-a261-21dd0e54cd20
md"
Oops! What happened here?!?

...Our custom feature can't actually encode values for tungsten or sulfur. How could we have known this ahead of time? There's a function for that!
"

# ╔═╡ d8025cba-e570-4da1-9f95-306b4be3b9a1
encodable_elements(fzn)

# ╔═╡ a1733b0d-859b-4621-b2a9-fc42a2afacc4
md"Let's try a different featurization; we'll replace `MeaningOfLife` with electronegativity."

# ╔═╡ 95114292-20c4-4973-a2fa-7064948bfced
new_fzn = GraphNodeFeaturization(["Block", "Atomic mass", "X"])

# ╔═╡ d982e350-c502-4c55-8043-bcefe10a180f
encodable_elements(new_fzn)

# ╔═╡ 4dc19543-6c97-4823-b2b5-559550654f6b
md"That looks better!"

# ╔═╡ b342a83e-740d-4b67-aafb-186827f03a6c
featurized_WS2 = featurize(WS2, new_fzn)

# ╔═╡ 33feaba7-dca9-4603-aa52-cd5bd3e7b4ed
propertynames(featurized_WS2)

# ╔═╡ ce3b6b9c-a588-48fd-b3b8-3c6b1a98ac22
md"As promised, we can still decode!"

# ╔═╡ e7bdb5d0-46ff-4244-89e4-0dbb518507a0
decode(featurized_WS2)

# ╔═╡ Cell order:
# ╟─8cc13d04-83a9-11ec-2d38-47d50cf23b83
# ╠═50fbffca-95d9-4b2b-955e-b5453574ddaf
# ╟─db3022ee-9241-4446-87ad-92c8359d1e31
# ╟─96095ce5-6edb-46f0-9884-397b318ed631
# ╠═0f71df51-9e39-4768-8965-d79c0105f3e8
# ╠═84d37c16-1afc-4ad9-9ae5-8beb5c8d3bc8
# ╟─b2263a81-f2f0-4b0a-a6f0-393e36c1e0ba
# ╠═8add5136-1a19-440c-a2d7-0d19f2d5e809
# ╠═fc212e05-bc3b-484c-84e4-e051f7c555df
# ╟─3fb29959-a7db-4e48-a527-e81ac85b486d
# ╠═718c6c24-05ff-4aca-bb54-787992e1a517
# ╟─a501a4db-fed5-43e7-bbc2-bf5657475ee1
# ╠═114a2eff-acef-42b1-bedc-96dfc260e127
# ╠═084f5ad9-6381-4663-97f6-f5f73ac2944d
# ╟─ee8549d0-6ba2-4749-be52-74814756deca
# ╟─6af1e757-79ce-4e5c-88f1-45b4d39a993f
# ╟─75b26589-b74f-4508-8d62-fbf40423d64a
# ╠═39e9b99a-3964-4f5d-a5a7-4bf3726c4429
# ╟─9ffc63ac-04cd-4524-88ad-2339035c8daf
# ╠═12e7461c-af53-4449-a800-4663fb4839ee
# ╠═75100b89-f2cc-40bc-9d8a-36f52831292e
# ╟─a9be0242-65ff-465d-8c1e-6e3792e38147
# ╟─0ee1995b-77a8-4415-b686-f876fa2761f8
# ╠═ee15cd74-bdc2-4437-a100-1d5857f33b7e
# ╟─4f73e95c-ffdd-4227-bbe5-14b80413cd46
# ╠═773c4a97-c5ab-486e-8286-08b3b7167e7d
# ╟─26eb23af-d1b2-43be-8f04-f820f441c8e3
# ╠═25909b37-9a83-44a3-a287-a19c7e542391
# ╟─22fe0263-ffc3-4cbe-9f7b-29f3a0ce6d97
# ╠═da01488b-efd0-4653-a9a2-b8defe93ff17
# ╟─bb6bcb31-2905-49af-8070-bb493069cc84
# ╠═644b1a5a-65e7-4731-8fd8-058670636b1d
# ╟─2ac9fd2a-44c6-4c40-b8dc-9ce47c8eb333
# ╠═1127b60a-b52d-4d30-8d4a-06a38dba4d04
# ╠═d6205d96-6cb3-438b-8adb-b5d6b8acaa45
# ╟─32fa1fd9-55cd-4615-a2bc-2a9e0b2524c1
# ╠═47cd8c6a-32ba-4f23-96bf-8b678c8b1721
# ╟─99a39408-d1f8-4afb-8f18-7bcd1b7faf08
# ╠═f9d92884-cdce-4b9c-984f-3f5bd69ea020
# ╠═8ed6f27e-d40f-4d05-a92c-3515980aad59
# ╟─bbcca9da-3027-4025-885f-21b12454a5aa
# ╟─b00a2447-06bc-4d93-b717-96f01771c176
# ╠═b94b89c1-5496-41d5-be02-6136763f1153
# ╠═b1595f7e-21c1-44d7-b942-2d652cb5a28b
# ╟─51db299f-5a1c-431c-a248-42be75a1e18f
# ╠═e436c12d-c5fe-4904-be9c-6802559a5fd1
# ╠═1093a1fb-edfb-4176-b310-4e1fa5572473
# ╟─6d75ebb9-3212-4058-8236-8c900b961019
# ╠═92c9954b-72e7-4d91-988a-915b2c5999d5
# ╠═91afd9c2-fc59-4b2b-b1f0-b9cf61a06ade
# ╠═2c849e35-9161-4ad4-aa80-6ee203205810
# ╟─c4a25ad0-ee4f-43f9-9040-5d1d39aaa72a
# ╠═317cd6d0-18b7-4fcd-8137-46d6ed2bdbff
# ╠═c307dc74-879e-45d7-8cd1-076580a66c43
# ╟─6e61b457-12dd-4038-8e4c-23bc3cc62ba5
# ╟─8ccfc9bc-0b4b-4ea7-89c7-86b7d3b756e4
# ╠═c10632e4-4f50-42a1-93e7-e4736ab0237e
# ╟─51d2f6d9-681e-46a4-89d8-d0d25a3fdebd
# ╟─165f0c92-0f01-4525-a678-51178f36c328
# ╠═498144b5-e599-45da-97b7-a65931429ab6
# ╟─e8d10af1-bde0-49ee-ac9a-33c82da7cd1c
# ╠═9faaa0b2-ef70-470f-9011-83a1367e6963
# ╟─44a3a9bb-7c6b-4270-bea8-71ecb3fa7aed
# ╠═96bb6f6d-ada0-4fb8-bdf2-f12de1e68b71
# ╟─2a619e89-8398-48e7-afb3-ba11a9557820
# ╠═06c294c4-fa86-4d17-891d-39af1c81586e
# ╟─d5535021-8c23-4640-a261-21dd0e54cd20
# ╠═d8025cba-e570-4da1-9f95-306b4be3b9a1
# ╟─a1733b0d-859b-4621-b2a9-fc42a2afacc4
# ╠═95114292-20c4-4973-a2fa-7064948bfced
# ╠═d982e350-c502-4c55-8043-bcefe10a180f
# ╟─4dc19543-6c97-4823-b2b5-559550654f6b
# ╠═b342a83e-740d-4b67-aafb-186827f03a6c
# ╠═33feaba7-dca9-4603-aa52-cd5bd3e7b4ed
# ╟─ce3b6b9c-a588-48fd-b3b8-3c6b1a98ac22
# ╠═e7bdb5d0-46ff-4244-89e4-0dbb518507a0

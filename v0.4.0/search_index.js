var documenterSearchIndex = {"docs":
[{"location":"types/codecs/#Codecs","page":"Codec","title":"Codecs","text":"","category":"section"},{"location":"types/codecs/","page":"Codec","title":"Codec","text":"Codec.OneHotOneCold","category":"page"},{"location":"types/codecs/#ChemistryFeaturization.Codec.OneHotOneCold","page":"Codec","title":"ChemistryFeaturization.Codec.OneHotOneCold","text":"OneHotOneCold(encode_f, decode_f, nbins, logspaced)\n\nAbstractCodec type which uses a dummy variable (as defined in statistical literature), i.e., which employs one-hot encoding and a one-cold decoding scheme.\n\n\n\n\n\n","category":"type"},{"location":"types/feature_descriptors/#Feature-Descriptors","page":"Feature Descriptors","title":"Feature Descriptors","text":"","category":"section"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"Feature descriptors store all necessary information to encode and decode feature values on various parts of an atoms object and appropriately combine them into a single object (vector, matrix, etc.) describing the value/values of the feature for the entire object.","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"For example, if an ElementFeatureDescriptor encodes a vector for each atom in an object, they could be concatenated together into a matrix with a column for each atom to describe a structure.","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"The type hierarchy of these objects is currently:","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"|---- AbstractFeatureDescriptor\n    |---- AbstractAtomFeatureDescriptor\n        |==== ElementFeatureDescriptor\n        |==== SpeciesFeatureDescriptor\n    |---- AbstractPairFeatureDescriptor\n        |==== PairFeatureDescriptor\n        |---- BondFeatureDescriptor\n            |==== BondType\n            |==== InRing\n            |==== IsConjugated","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"where  ---- = Abstract Type and ==== = Concrete Type","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"More details on each of these types is below, and more types (e.g. environment features) will be implemented in the future!","category":"page"},{"location":"types/feature_descriptors/#Functionality-common-to-all-feature-descriptors","page":"Feature Descriptors","title":"Functionality common to all feature descriptors","text":"","category":"section"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"they should be callable on atoms objects and return encoded features\nSimilarly, decode should work...","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"decode(::FeatureDescriptor.AbstractFeatureDescriptor, ::Any)","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"the function encodable_elements should be defined on all feature descriptors, as it will be used to verify that a feature can be encoded for every atom in a structure.","category":"page"},{"location":"types/feature_descriptors/#Atom-Feature-Descriptors","page":"Feature Descriptors","title":"Atom Feature Descriptors","text":"","category":"section"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"These types encode features for single atoms in a structure. The abstract parent type is AtomFeatureDescriptor.","category":"page"},{"location":"types/feature_descriptors/#Element-Feature-Descriptors","page":"Feature Descriptors","title":"Element Feature Descriptors","text":"","category":"section"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"An ElementFeatureDescriptor's encoded values are defined only by the elemental identity of an atom. Examples include atomic mass and block (s, p, d, or f) in the periodic table.","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"ElementFeatureDescriptor","category":"page"},{"location":"types/feature_descriptors/#ChemistryFeaturization.FeatureDescriptor.ElementFeatureDescriptor","page":"Feature Descriptors","title":"ChemistryFeaturization.FeatureDescriptor.ElementFeatureDescriptor","text":"ElementFeatureDescriptor\n\nDescribe features associated with individual atoms that depend only upon their elemental identity\n\nFields\n\nname::String: Name of the feature\nencoder_decoder::AbstractCodec: Codec defined which handles the feature's encoding and decoding logic\ncategorical::Bool: flag for whether the feature is categorical or continuous-valued\nlookup_table::DataFrame: table containing values of feature for every encodable element\n\n\n\n\n\n","category":"type"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"In the example, below, we encode the block of each atom in a hydrogen molecule. The result is two hcatted vectors [1 0 0 0], indicating hydrogen is s-block.","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"TODO: check that this test passes once new version is tagged","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"H2 = AtomGraph([0. 1.; 1. 0.], [\"H\", \"H\"])\nblock = ElementFeatureDescriptor(\"Block\")\nblock(H2)\n\n# output\n4×2 Matrix{Float64}:\n 1.0  1.0\n 0.0  0.0\n 0.0  0.0\n 0.0  0.0","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"Because they are defined only by the element, values for these features can be tabulated in a lookup table. Many commonly-desired element features are included in the atom_data_df DataFrame, but you can also define custom lookup tables for other features by utilizing the lookup_table keyword of the ElementFeatureDescriptor constructor.","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"TODO: add remark about encoding options once that PR is merged","category":"page"},{"location":"types/feature_descriptors/#Species-Feature-Descriptor","page":"Feature Descriptors","title":"Species Feature Descriptor","text":"","category":"section"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"A SpeciesFeatureDescriptor's encoded values depend on its local environment. Examples are an atom's format oxidation state, or whether it is part of an aromatic ring.","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"TODO: more details once we have better examples","category":"page"},{"location":"types/feature_descriptors/#Pair-Feature-Descriptors","page":"Feature Descriptors","title":"Pair Feature Descriptors","text":"","category":"section"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"Pair feature descriptors encode features of pairs of atoms. The abstract parent type is AbstractPairFeatureDescriptor.","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"The concrete type PairFeatureDescriptor encodes information about any pair of atoms, such as the distance between them.","category":"page"},{"location":"types/feature_descriptors/#Bond-Feature-Descriptors","page":"Feature Descriptors","title":"Bond Feature Descriptors","text":"","category":"section"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"Bond feature descriptors are defined only for two atoms that are bonded to each other.","category":"page"},{"location":"types/feature_descriptors/","page":"Feature Descriptors","title":"Feature Descriptors","text":"TODO: more details here","category":"page"},{"location":"types/featurizations/#Featurization","page":"Featurization","title":"Featurization","text":"","category":"section"},{"location":"contributing/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"We welcome community contributions! This page contains some guidance to make the process go smoothly.","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"Pages = [\"contributing.md\"]","category":"page"},{"location":"contributing/#General-Guidelines","page":"Contributing","title":"General Guidelines","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"For some high-level general guidance, please see CONTRIBUTING.md on the GitHub repo.","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"The remainder of this page includes some pointers about package structure, etc. that could be helpful depending on what sorts of functionality you're interested in adding to the package. An understanding of \"what goes where\" will help in making sure you put your code in the right place! We also strongly encourage you to read the Terminology/Philosophy page, as it will help to understand the best ways to add different types of functionality.","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"TODO: flesh out everything below","category":"page"},{"location":"contributing/#Implementing-new-atoms-objects","page":"Contributing","title":"Implementing new atoms objects","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"does it need to be a new type? Could it be accommodated by an existing type, or by an existing type with an expansion of functionality?\nmake sure to subtype AbstractAtoms\ncode belongs in src/atoms/, export statements belong in...\nmake sure any/all sensible FD's and fzn's dispatch on it appropriately","category":"page"},{"location":"contributing/#Implementing-new-feature-descriptors","page":"Contributing","title":"Implementing new feature descriptors","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"does it need to be a new type? Could it be accommodated by an existing type, or by an existing type with an expansion of functionality?\nmake sure to place it in the right place in the type hierarchy\ncode belongs in src/features/\nmake sure it encodes on as many atoms objects as are sensible","category":"page"},{"location":"contributing/#Implementing-new-featurization-schemes","page":"Contributing","title":"Implementing new featurization schemes","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"does it need to be a new type? Could it be accommodated by an existing type, or by an existing type with an expansion of functionality?\nshould subtype AbstractFeaturization\ncode belongs in src/featurizations/\nmake sure it encodes on as many atoms objects as are sensible","category":"page"},{"location":"changelog/#Changelog","page":"Changelog","title":"Changelog","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"I'm generally trying to adhere to semver here. This means that in v0.*, assume breaking changes are always possible, even without a major version bump...however, I will try to always note them here if they happen...","category":"page"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"Categories to include for each release, if relevant: breaking, added, fixed, removed/deprecated","category":"page"},{"location":"changelog/#v0.3.0","page":"Changelog","title":"v0.3.0","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"Major restructure to make extensibility far easier and allow better sharing of elements between featurization schemes! Also streamlined the codebase a lot by organizing things into sensible modules, so technically basically everything is a breaking change because all of the imports are from the submodules now.","category":"page"},{"location":"changelog/#Breaking/Added","page":"Changelog","title":"Breaking/Added","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"create modules for different functionalities: Atoms, Features, Featurization, Utils, etc.\nseparate out encoding/decoding functionality into a new AbstractCodec type and associated Codec module\nmerged functionality of AtomGraph and WeaveMol into AtomGraph object, with a more generic encoded_features field","category":"page"},{"location":"changelog/#v0.2.2-[2021-02-22]","page":"Changelog","title":"v0.2.2 [2021-02-22]","text":"","category":"section"},{"location":"changelog/#Added","page":"Changelog","title":"Added","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"updated pretty printing for AtomGraph to include id field","category":"page"},{"location":"changelog/#Fixed","page":"Changelog","title":"Fixed","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"added some docstrings, updated others","category":"page"},{"location":"changelog/#v0.2.1-[2021-02-22]","page":"Changelog","title":"v0.2.1 [2021-02-22]","text":"","category":"section"},{"location":"changelog/#Breaking","page":"Changelog","title":"Breaking","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"rename build_atom_feats to build_featurization to avoid ambiguity with make_feature_vectors","category":"page"},{"location":"changelog/#Fixed-2","page":"Changelog","title":"Fixed","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"proper docstrings for more things incl. AtomFeat, AtomGraph, add_features!","category":"page"},{"location":"changelog/#v0.2.0-[2021-02-16]","page":"Changelog","title":"v0.2.0 [2021-02-16]","text":"","category":"section"},{"location":"changelog/#Breaking-2","page":"Changelog","title":"Breaking","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"add id field to AtomGraph (should only be breaking for reading in serialized graphs)\nbuild_graphs_batch now returns list of graphs and optionally serializes rather than always serializing and never returning","category":"page"},{"location":"changelog/#Added-2","page":"Changelog","title":"Added","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"add_features_batch! function\nAtomGraph constructor directly from adjacency matrix (previously was only from a SimpleWeightedGraph)","category":"page"},{"location":"changelog/#Fixed-3","page":"Changelog","title":"Fixed","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"remove deprecated syntax of opening DataFrames from CSV\nmade separate files into modules to avoid redefinition warnings during precompilation","category":"page"},{"location":"changelog/#v0.1.1-[2021-02-10]","page":"Changelog","title":"v0.1.1 [2021-02-10]","text":"","category":"section"},{"location":"changelog/#Fixed-4","page":"Changelog","title":"Fixed","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"add check for NaN values in graph laplacian","category":"page"},{"location":"changelog/#v0.1.0-[2020-12-22]","page":"Changelog","title":"v0.1.0 [2020-12-22]","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"Initial release!","category":"page"},{"location":"changelog/#Added-3","page":"Changelog","title":"Added","text":"","category":"section"},{"location":"changelog/","page":"Changelog","title":"Changelog","text":"Create AtomGraph and AtomFeat types\nBasic graph visualization functions\nGraph-building from CIF files via the cgcnn.py \"cutoff\" method, with support for nonperiodic systems as well","category":"page"},{"location":"types/atoms/#Atoms-Objects","page":"Atoms Objects","title":"Atoms Objects","text":"","category":"section"},{"location":"types/atoms/","page":"Atoms Objects","title":"Atoms Objects","text":"Atoms objects (terminology borrowed from ASE) store information about the structure of a molecule, crystal, etc. as well as, optionally, encoded features and the featurization used to encode them. The parent abstract type is AbstractAtoms.","category":"page"},{"location":"types/atoms/#AtomGraph","page":"Atoms Objects","title":"AtomGraph","text":"","category":"section"},{"location":"types/atoms/","page":"Atoms Objects","title":"Atoms Objects","text":"The AtomGraph type is used to store atomic graph representations.","category":"page"},{"location":"types/atoms/","page":"Atoms Objects","title":"Atoms Objects","text":"Atoms.AtomGraph","category":"page"},{"location":"types/atoms/#ChemistryFeaturization.Atoms.AtomGraph","page":"Atoms Objects","title":"ChemistryFeaturization.Atoms.AtomGraph","text":"AtomGraph\n\nA type representing an atomic structure as a graph (gr).\n\nFields\n\ngraph::SimpleWeightedGraph{<:Integer,<:Real}: the graph representing the structure. See build_graph for more on generating the weights.\nelements::Vector{String}: list of elemental symbols corresponding to each node of the graph\nlaplacian::Matrix{<:Real}: Normalized graph Laplacian matrix, stored to speed up convolution operations by avoiding recomputing it every pass.\nid::String: Optional, an identifier, e.g. to correspond with tags/labels of an imported dataset.\n\n\n\n\n\n","category":"type"},{"location":"types/atoms/#WeaveMol","page":"Atoms Objects","title":"WeaveMol","text":"","category":"section"},{"location":"#ChemistryFeaturization.jl","page":"Home","title":"ChemistryFeaturization.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Data types and featurization schemes for Chemellia models and beyooond!","category":"page"},{"location":"","page":"Home","title":"Home","text":"ChemistryFeaturization.jl is meant to be a unified interface for translating atomic structures (molecules, crystals, etc.) into data structures and sets of features to be used in machine learning models provided by the Alchemy suite of packages, such as AtomicGraphNets.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package is in development as part of the ACED project, funded by ARPA-E DIFFERENTIATE and coordinated by a team from Carnegie Mellon University, in collaboration with Julia Computing, Citrine Informatics, and MIT. Dr. Rachel Kurchin is the lead developer.","category":"page"},{"location":"#Purpose","page":"Home","title":"Purpose","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is intended to serve as a \"helper package\" of sorts for other packages that actually build the models, and NOT as a standalone tool (for instance, ChemistryFeaturization.jl it does not itself implement any models).","category":"page"},{"location":"","page":"Home","title":"Home","text":"It provides flexible, modular data types and functions for starting from a \"bare\" atomic structure of a molecule or crystal, featurizing it with various tabulated or computed properties of its constituent atoms, bonds, etc., and encoding those features into a format appropriate to ingest into a machine learning model. Critically, it also provides functions for inverting this process, i.e. decoding the featurization, as all necessary metadata is retained.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: CF_flowchart)","category":"page"},{"location":"","page":"Home","title":"Home","text":"To learn more about the package, it would be good to start with the Terminology/Philosophy section to get a handle on the vocabulary as we use it, and also some insights about the design philosophy behind the package. Then go ahead and browse the sidebar for more!","category":"page"},{"location":"terminology/#Terminology/Philosophy","page":"Terminology","title":"Terminology/Philosophy","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"There are a lot of seemingly similar terms used for quantities in this package that refer to disparate things (or, are used slightly differently by other people in other places). Here, we try to best define these terms as we intend them. Further down, once the terms are defined, we elaborate on why the package is designed the way it is.","category":"page"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"Pages = [\"terminology.md\"]\nDepth = 3","category":"page"},{"location":"terminology/#General-Terms","page":"Terminology","title":"General Terms","text":"","category":"section"},{"location":"terminology/#Feature","page":"Terminology","title":"Feature","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"A quality or quantity associated with an atom that we wish to encode, such as atomic mass, row in the periodic table, etc.","category":"page"},{"location":"terminology/#Encoding","page":"Terminology","title":"Encoding","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"The process of translating the value of a feature from its human-readable form (such as a float or a string) to whatever form will be ingested by a machine learning model. This could be as simple as an equality operation, but more often is, e.g. building a one-hot vector.","category":"page"},{"location":"terminology/#Decoding","page":"Terminology","title":"Decoding","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"The inverse process to encoding. Note that in many cases (e.g. a continuous-valued feature encoded to a one-hot vector), the process isn't fully invertible, i.e. you can't get back a precise value but rather only a range corresponding to the associated onehot bin.","category":"page"},{"location":"terminology/#Data-types-in-ChemistryFeaturization","page":"Terminology","title":"Data types in ChemistryFeaturization","text":"","category":"section"},{"location":"terminology/#Feature-Descriptor","page":"Terminology","title":"Feature Descriptor","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"Describes the \"features of a feature\" – i.e. its name, possible values, instructions for encoding it, etc., but does NOT store an actual instance of its value.","category":"page"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"For more on the available types of feature descriptors, see Feature Descriptors.","category":"page"},{"location":"terminology/#AbstractCodec","page":"Terminology","title":"AbstractCodec","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"Component of a feature descriptor that stores the actual encoding/decoding functions. ","category":"page"},{"location":"terminology/#Atoms-Object","page":"Terminology","title":"Atoms Object","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"Describes a molecule, crystal, etc. in whatever representation will be ingested by an ML model (e.g. a graph), and can also store encoded features of that structure. ","category":"page"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"For more on the available types of atoms objects, see Atoms Objects.","category":"page"},{"location":"terminology/#Featurization-Object","page":"Terminology","title":"Featurization Object","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"Stores sets of feature descriptors and instructions for combining the values they encode on an atoms object.","category":"page"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"For more on the available types of featurization objects, see Featurization.","category":"page"},{"location":"terminology/#Design-Philosophy","page":"Terminology","title":"Design Philosophy","text":"","category":"section"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"Points to elaborate on here...","category":"page"},{"location":"terminology/","page":"Terminology","title":"Terminology","text":"maintaining transparency/decodability, as well as user's choice to include (or not) particular features, and choose how they are encoded, as opposed to a \"black-box\" scheme with little to no customizability\nseparation of concerns/modularity...as many things should be \"plug-and-play\" with each other as possible (e.g. swapping in different codecs to FD's, different FD's to featurizations, different featurizations to atoms objects)\nFD's and atoms objects are generic and fairly reusable across models, featurizations are less so (closer to one-to-one relationship between a featurization type and a model type, but still have flexibility to easily include/exclude different features)","category":"page"}]
}

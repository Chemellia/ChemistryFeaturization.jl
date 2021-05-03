module ChemistryFeaturization

using SimpleWeightedGraphs
using Reexport

# define all the abstract types 
export AbstractAtoms, AbstractFeature, AbstractFeaturization
abstract type AbstractAtoms end
abstract type AbstractFeature end
abstract type AbstractFeaturization end

include("utils/Utils.jl")
@reexport using .Utils.AtomFeatureUtils

#= ATOMS OBJECTS

Should have a field called `elements` that is a list of strings of elemental symbols.

AtomGraph shouldn't have to change too much from its current incarnation. WeaveMol is 
not currently a thing that can be directly fed into a model, so that will have to get
updated. Both will need some way to specify which types of features can be attached 
to them, maybe this is a place for the holy trait design pattern?

=#

# link to guidance in docs about how to implement new feature types

# export...
export AtomGraph, visualize
#export WeaveMol#, ...

# include...
include("atoms/atomgraph.jl")
#include("weavemol.jl")


#= FEATURE OBJECTS

I decided on an abstract type to allow the dispatch shown below in the "magical encoding"
and "magical decoding" bits.

Representative examples – feature: feature object type, input to encoder, output of encoder:

block (s,p,d,f):            AtomFeature (contextual=false), String, Flux.OneHotMatrix 
electronegativity (binned): AtomFeature (contextual=false), Float32, Flux.OneHotMatrix
electronegativity (direct): AtomFeature (contextual=false), Float32, Vector{Float32}
oxidation state:            AtomFeature (contextual=true) Int, Vector{Float32}
distance between atoms:     PairFeat, Float32, Matrix{Float32}
bond type:                  PairFeat, String, Array{Float32,3}

(bond type is categorical and gets one-hot encoded, so the first two indices of the Matrix
should be atom indices and the third should be indexing into the one-hot vector, which should just contain all zeros if the two atoms are not bonded (or we add a bin to the one-hot encoding to indicate that)

All subtypes should define `encode_f` and `decode_f`
=#

# link to guidance in docs about how to implement new feature types

# export...
export AtomFeature, PairFeature

# include...
include("features/atomfeature.jl")
include("features/pairfeature.jl")

# generic encode
# docstring
function (f::AbstractFeature)(a::AbstractAtoms)
    f.encode_f(a)
end

export decode
# generic decode
# docstring
function decode(f::AbstractFeature, encoded_feature)
    f.decode_f(encoded_feature)
end

#= FEATURIZATION OBJECTS
All such objects should define at least one list of <:AbstractFeature objects and either work according to the generic featurize! defined herein or dispatch featurize! if customized behavior is needed.
=#

# export...
export GraphNodeFeaturization, WeaveFeaturization, featurize!

# include...
include("featurizations/graphnodefeaturization.jl")
include("featurizations/weavefeaturization.jl")

# generic featurize!
# this assumes that `a` has fields with names corresponding to each field in `fzn`, if not you need to dispatch this function to your specific case
# TODO: maybe add option to exclude field names from iteration over fzn?
# docstring
function featurize!(a::AbstractAtoms, fzn::AbstractFeaturization)
    # loop over fields in featurization, each one is a list of features
    # encode each feature in that list and assign the results to the
    # field of the same name in `a`
    for feats_list in fieldnames(fzn)
        encoded = reduce(vcat, map((x)->x(a), feats_list))
        setproperty!(a, feats_list, encoded)
    end
    a.featurization = fzn
end


# NEXT: 
# building features
# building featurizations
# featurize atomgraphs
# update tests


# old stuff for archival purposes for now...

# export AtomFeature, atom_data_df, build_featurization, make_feature_vectors, decode_feature_vector, default_nbins
# include("atomfeat.jl")

# export AtomGraph, normalized_laplacian, add_features!, add_features_batch!, visualize_graph
# include("atomgraph.jl")

# export inverse_square, exp_decay, build_graph, build_graphs_batch, read_graphs_batch
# include("pmg_graphs.jl")
# using .graph_building: inverse_square, exp_decay, build_graph, build_graphs_batch, read_graphs_batch

# # TODO: possibly move all the Weave stuff to another package altogether, if not need to tidy up modules/exports
# export smiles_atom_features, smiles_bond_features
# include("weave_fcns.jl")
# using .weave_fcns: smiles_atom_features, smiles_bond_features, chem

# export weave_featurize
# include("featurize.jl")
# using .featurize: weave_featurize

end

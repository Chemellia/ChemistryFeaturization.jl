#= 
Vocabulary for the purposes of this document (I am VERY open to clearer/more 
precise terminology suggestions as I'm not necessarily happy with these but 
rather needed to set some conventions of what means what):

FEATURE OBJECT: a struct that describes "the features of a feature" – i.e. its name,
                possible values, how to encode it, etc. but does NOT store an actual
                instance of its value

ATOMS OBJECT: (borrowed from ASE terminology) a struct that describes a molecule,
              crystal, etc. in whatever representation will be ingested by an ML
              model (e.g. a graph), and can also store encoded features. If it stores
              encoded features, it MUST include requisite metadata (typically in the
              form of a list of FEATURE OBJECTs and potentially also a function handle)
              to decode the function again. Examples incldue AtomGraph, WeaveMol (though
              their current incarnations may need to be adapted)

ENCODING: the process of translating a feature from its human-readable form (e.g. a 
          float, string, etc.) to whatever will be ingested by ML model (could be as
          simple as making a copy, more often is, e.g. building a one-hot vector)

DECODING: The inverse process to encoding. Note that in many cases (e.g. a continuous-
          valued feature encoded to a one-hot vector), the process isn't fully invertible,
          i.e. you can't get back a precise value but rather only a range corresponding
          to the associated onehot bin

FEATURIZATION: the process of attaching encoded features to an atoms object (in a 
               slight abuse of terminology, a featurization function  may also perform
               feature encoding).

File organization: probably eventually will be cleanest to have a `feat` folder with 
separate files defining each feature type and similarly for the atoms object types and
featurizations.

In addition, there should be some place to drop functions for automatically building 
features and featurizations for commonly desired features (i.e. things from lookup table
and things making use of that dictionary of rdkit functions).

For now, I'm sketching this all out in one file though.
=#

#= FEATURE OBJECTS

I decided on an abstract type to allow the dispatch shown below in the "magical encoding"
and "magical decoding" bits.

Tn is the "natural" type of feature values, Te is the type of encoded values for the WHOLE
STRUCTURE (not e.g. for a single atom).

Representative examples – feature: feature object type, Tn, Te:

block (s,p,d,f):            AtomFeat (contextual=false), String, Flux.OneHotMatrix 
electronegativity (binned): AtomFeat (contextual=false), Float32, Flux.OneHotMatrix
electronegativity (direct): AtomFeat (contextual=false), Float32, Vector{Float32}
oxidation state:            AtomFeat (contextual=true) Int, Vector{Float32}
distance between atoms:     PairFeat, Float32, Matrix{Float32}
bond type:                  PairFeat, String, Array{Float32,3}

(bond type is categorical and gets one-hot encoded, so the first two indices of the Matrix
should be atom indices and the third should be indexing into the one-hot vector, which should just contain all zeros if the two atoms are not bonded (or we add a bin to the one-hot encoding to indicate that)

All subtypes should define `encode_f` and `decode_f`
`encode_f` should take in an atoms object and return Te
`decode_f` should take in something of type Te and return something of type Tn
=#
abstract type AbstractFeature{Tn,Te} end

# I think we should be able to do this too...there might be a more idiomatic way
# to do what I'm showing below, perhaps via an "actual" holy trait type thing
isatoms(::Type{AtomGraph}) = true
isatoms(::Type{WeaveMol}) = true
# add more for other types later on...
isatoms(::Type{Any}) = false

# magical encoding
function (f<:AbstractFeature{Tn,Te})(atomsobj::T) where T
    @assert isatoms(T) "Features can only be computed on atomic structures!"
    f.encode_f(atomsobj)
end

# magical decoding
(f<:AbstractFeature{Tn,Te})(encoded::Te) where {Tn,Te} = f.decode_f(encoded)

#=
Feature of a single atom.

(logspaced field not necessary because it can be embedded in encode_f and decode_f through 
keywords to constructor (see below))
=#
struct AtomFeat{Tn,Te}<:AbstractFeature{Tn,Te}
    name # symbol (for easy DataFrame indexing) or String (for easy other things)?
    encode_f
    decode_f
    categorical::Bool
    contextual::Bool # can get from elemental lookup table (false) or not (true)?
    length::Int # length of encoded vector, useful for computing eventual sizes of things probably?
end

#=
we'll define a bunch of automatic stuff for building AtomFeats with built-in data
that will essentially copy the current contents of atomfeat.jl, ideally with some
additions enabling the user to augment the lookup table

should also have convenience functions for building encode_f and decode_f via keywords for, e.g.:
    choosing one-hot (and how many bins, etc.) vs. direct float encoding
    logspaced vs. linear spaced binning for one-hot encoding
    maybe other stuff too
...
=#

# then basically copy a bunch of the weave stuff, plus potentially things like oxidation state etc. from pymatgen belong here too

#=
Feature of a pair of atoms. Currently only used in WeaveModel, but may eventually
have a version of AtomGraph that allows edge features...
=#
struct PairFeat{Tn,Te}<:AbstractFeature{Tn,Te}
    name
    encode_f
    decode_f
    length::Int # maybe, maybe not (does constrain/assume vector Te)
    # probably needs some other stuff...
end
# same idea here...copy from weave stuff...
# a note: I think it's more elegant (and probably efficient?) if pair features assert
# always arrays where the first two indices are atom indices and any further are indexing
# into the feature if it is itself higher-dimensional. This is not how the weave featurizer
# currently works...

#=
This is just an example of another sort of feature one might add, but it's not needed/used
currently.
=#
struct StructureFeat{Tn,Te}<:AbstractFeature{Tn,Te}
    name
    encode_f
    decode_f
    # etc...
end

#= FEATURIZATION OBJECTS
All such objects should define at least one list of <:AbstractFeature objects and a `combine`
function that specifies how to put them together for eventual attachment to an atoms 
object. 
=#

abstract type AbstractFeaturization end

struct GraphNodeFeaturization <: AbstractFeaturization
    atom_feats::Vector{AtomFeat}
    feature_vectors::Dict{String,Vector{Float32}} # map from element symbol to vector
    combine
end

struct WeaveFeaturization <: AbstractFeaturization
    atom_feats::Vector{AtomFeat}
    pair_feats::Vector{PairFeat}
    combine
end


#= ATOMS OBJECTS

AtomGraph shouldn't have to change too much from its current incarnation. WeaveMol is 
not currently a thing that can be directly fed into a model, so that will have to get
updated. Both will need some way to specify which types of features can be attached 
to them, maybe this is a place for the holy trait design pattern?

Example: AtomGraph can take ElementFeat and ComputedAtomFeat but not PairFeat or StructureFeat, WeaveMol can take ElementFeat, ComputedAtomFeat, and PairFeat but also not StructureFeat
=#

abstract type AbstractAtoms end

# it is kind of cool to subtype this from LightGraphs so all that machinery "just works"
# that being said, all of it would "just work" by just pulling out ag.graph so may be
# more cool than actually practical, hence I've done the abstract type here
mutable struct AtomGraph <: AbstractAtoms
    graph::SimpleWeightedGraph{Int32,Float32}
    elements::Vector{String}
    lapl::LightGraphs.LinAlg.NormalizedLaplacian
    features::Matrix{Float32} # if we add edge features this type will have to relax
    featurization::GraphNodeFeaturization
    id::String # or maybe we let it be a number too?
end

mutable struct WeaveMol <: AbstractAtoms
    smiles::String
    elements::Vector{String}
    features::Tuple{SomethingOrOther} # I need to look more carefully to figure this out heh
    featurization::WeaveFeaturization
    id # probably makes sense to have this in addition to smiles? Maybe?
end

# constructors etc. as before
# maybe some cutesy stuff like dispatching things like length as length(elements)
# pretty printing stuff

# generic featurize...there's probably a better way to write this...
function featurize!(a<:AbstractAtoms, f<:AbstractFeaturization)
    feats_lists = [field for field in fieldnames(typeof(f)) if field!==:combine]
    encoded_features = Dict(feats_list => [] for feats_list in feats_lists)
    for feats_list in feats_lists
        features = map(feat->feat(a), feats_list)
        encoded_features[feats_list] = features
    end
    a.features = f.combine(encoded_features)
    a.featurization = f
end

# then there would be specific versions of featurize! for particular featurizations
# in order to facilitate other function signatures, etc.

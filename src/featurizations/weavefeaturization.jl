# some commentary

using ..ChemistryFeaturization.AbstractType: AbstractFeaturization
using ..ChemistryFeaturization.FeatureDescriptor:
    AbstractAtomFeatureDescriptor, AbstractPairFeatureDescriptor
using ..ChemistryFeaturization.Codec: DirectCodec
using ..ChemistryFeaturization.Utils

struct WeaveFeaturization <: AbstractFeaturization
    atom_features::Vector{<:AbstractAtomFeatureDescriptor}
    bond_features::Vector{<:BondFeatureDescriptor}
    pair_features::Vector{<:AbstractPairFeatureDescriptor}
end

function WeaveFeaturization(default_atom_feature_list = ["symbol",
                                                         "degree",
                                                         "implicithconnected",
                                                         "charge",
                                                         "radical_electrons",
                                                         "hybridization",
                                                         "isaromatic",
                                                         "hydrogenconnected"],
                            default_bond_feature_list = ["bondorder", "isaromaticbond", "isringbond"])
  species_feats = intersect(keys(Utils.SpeciesFeatureUtils.sfd_names_props), default_atom_feature_list)
  
  atoms = SpeciesFeatureDescriptor.(species_feats)
  bonds = BondFeatureDescriptor.(default_bond_feature_list)
  # pairs = PairFeatureDescriptor.(default_atom_feature_list)
  WeaveFeaturization(atoms, bonds, bonds)
end

function encodable_elements(fzn::WeaveFeaturization)
    intersect([encodable_elements(f) for f in fzn.atom_features]...),
    intersect([encodable_elements(f) for f in fzn.bond_features]...),
    intersect([encodable_elements(f) for f in fzn.pair_features]...)
end

# function encode(fzn::GraphNodeFeaturization, ag::AtomGraph)
#     encoded = reduce(vcat, map((x) -> encode(x, ag), fzn.features))
#     return encoded
# end

function atom_features(feat, mol; kw...)
  
end

const DEEPCHEM_ATOM_SYMBOLS = [
        "C",
        "N",
        "O",
        "S",
        "F",
        "Si",
        "P",
        "Cl",
        "Br",
        "Mg",
        "Na",
        "Ca",
        "Fe",
        "As",
        "Al",
        "I",
        "B",
        "V",
        "K",
        "Tl",
        "Yb",
        "Sb",
        "Sn",
        "Ag",
        "Pd",
        "Co",
        "Se",
        "Ti",
        "Zn",
        "H",  # H?
        "Li",
        "Ge",
        "Cu",
        "Au",
        "Ni",
        "Cd",
        "In",
        "Mn",
        "Zr",
        "Cr",
        "Pt",
        "Hg",
        "Pb",
        "Unknown"
      ]

# weave_mols = weave_featurize(smiles[Int((b_i-1)*batch_size+1):Int(b_i*batch_size)],
#                              atom_feature_list = atom_feature_list,
#                              bond_feature_list = bond_feature_list)
# y_batch = [Y[Int((b_i-1)*batch_size+1):Int(b_i*batch_size)]]
# default_atom_feature_list = ["symbol","degree","implicit_valence","formal_charge","radical_electrons","hybridization","aromaticity","total_H_num" ]
# default_bond_feature_list = ["bond_type","isConjugated","isInring"]

struct FeaturizedWeave
  atom_features
  bond_features
  pair_features
end

function encode(fzn::WeaveFeaturization, ag::AtomGraph; atom_feature_kwargs = (;),
                                                       bond_feature_kwargs = (;),
                                                       pair_feature_kwargs = (;))
  af = map(x -> encode(x, ag, atom_feature_kwargs...), fzn.atom_features) 
  bf = map(x -> encode(x, ag, bond_feature_kwargs...), fzn.bond_features)
  pf = map(x -> encode(x, ag, pair_feature_kwargs...), fzn.pair_features)
  FeaturizedWeave(af, bf, pf)
end

# function decode(fzn::GraphNodeFeaturization, encoded::Matrix{<:Real})
#     num_atoms = size(encoded, 2)
#     nbins = [output_shape(f) for f in fzn.features]
#     local decoded = Dict{Integer,Dict{String,Any}}()
#     for i = 1:num_atoms
#         #println("atom $i")
#         chunks = chunk_vec(encoded[:, i], nbins)
#         decoded[i] = Dict{String,Any}()
#         for (chunk, f) in zip(chunks, fzn.features)
#             #println("    $(f.name): $(decode(f, chunk))")
#             decoded[i][f.name] = decode(f, chunk)
#         end
#     end
#     return decoded
# end

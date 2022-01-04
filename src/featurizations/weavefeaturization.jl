# some commentary

using ..ChemistryFeaturization.AbstractType: AbstractFeaturization
using ..ChemistryFeaturization.FeatureDescriptor:
    AbstractAtomFeatureDescriptor, AbstractPairFeatureDescriptor
using ..ChemistryFeaturization.Codec: DirectCodec
using ..ChemistryFeaturization.Utils

struct WeaveFeaturization <: AbstractFeaturization
    element_features::Vector{<:AbstractAtomFeatureDescriptor}
    atom_features::Vector{<:AbstractAtomFeatureDescriptor}
    bond_features::Vector{<:BondFeatureDescriptor}
    pair_features::Vector{<:AbstractPairFeatureDescriptor}
end

function WeaveFeaturization(element_feature_list = ["Atomic no"],
                            species_feature_list = ["degree",
                                                    "implicithconnected",
                                                    "charge",
                                                    "radical_electrons",
                                                    "hybridization",
                                                    "isaromatic",
                                                    "hydrogenconnected"],
                            bond_feature_list = ["bondorder", "isaromaticbond", "isringbond"])
  elements = ElementFeatureDescriptor.(element_feature_list)  
  species = SpeciesFeatureDescriptor.(species_feature_list)
  bonds = BondFeatureDescriptor.(bond_feature_list)
  # pairs = PairFeatureDescriptor.(default_atom_feature_list)
  WeaveFeaturization(elements, species, bonds, bonds)
end

WeaveFeaturization(; kw...) = WeaveFeaturization(values(kw)...)

function encodable_elements(fzn::WeaveFeaturization)
    intersect([encodable_elements(f) for f in fzn.atom_features]...),
    intersect([encodable_elements(f) for f in fzn.bond_features]...),
    intersect([encodable_elements(f) for f in fzn.pair_features]...)
end

function atom_features(feat, mol; kw...)
  
end

function bond_feature(bond, mol; kw...)
  
end

function pair_features(feat, mol; kw...)

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
  sf = mapreduce(x -> encode(x, ag, atom_feature_kwargs...), vcat, fzn.atom_features)
  ef = mapreduce(x -> encode(x, ag, atom_feature_kwargs...), vcat, fzn.element_features)
  atom_and_elements = vcat(sf, ef)
  bf = cat(map(x -> encode(x, ag, bond_feature_kwargs...), fzn.bond_features)..., dims = 3)
  pf = cat(map(x -> encode(x, ag, pair_feature_kwargs...), fzn.pair_features)..., dims = 3)
  # Return FeaturizedAtoms here
  atom_and_elements, vcat(bf, pf)
  # FeaturizedWeave(atom_and_elements, bf, pf)
end

function Base.show(io::IO, fzn::WeaveFeaturization)
  println(io, "WeaveFeaturization(")
  println(io, "  Species Features: $(map(x -> x.name, fzn.element_features))")
  println(io, "  Atom Features: $(map(x -> x.name, fzn.atom_features))")
  println(io, "  Bond Features: $(map(x -> x.name, fzn.bond_features))")
  println(io, ")")
end

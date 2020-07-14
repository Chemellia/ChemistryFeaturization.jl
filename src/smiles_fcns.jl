# Sean to populate with ancillary functions for handling SMILES strings
using DataStructures
using Flux: onehot
using PyCall
chem = pyimport("rdkit.Chem")

# currently supported atomic features and bond features
const atom_feat_list = ["symbol","degree","implicit_valence","formal_charge","radical_electrons","hybridization","aromaticity","total_H_num" ]
const bond_feat_list = ["bond_type","isConjugated","isInring"]

# Map atomic feature names to RDKit functions
const rdkit_atom_function = Dict(
    "symbol"=>chem.Atom.GetSymbol,
    "degree"=>chem.Atom.GetDegree,
    "implicit_valence"=>chem.Atom.GetImplicitValence,
    "formal_charge"=>chem.Atom.GetFormalCharge,
    "radical_electrons"=>chem.Atom.GetNumRadicalElectrons,
    "hybridization"=>chem.Atom.GetHybridization,
    "aromaticity"=>chem.Atom.GetIsAromatic,
    "total_H_num"=>chem.Atom.GetTotalNumHs
)

# Map bond feature names to RDKit functions
const rdkit_bond_function = Dict(
    "bond_type"=>chem.Bond.GetBondType,
    "isConjugated"=>chem.Bond.GetIsConjugated,
    "isInring"=>chem.Bond.IsInRing
)

# bins for atomic features
const atom_feat_bins = Dict(
    "symbol" => [:C,:N,:O,:S,:F,:Si,:P,:Cl,:Br,:Mg,:Na,:Ca,:Fe,:As,:Al,:I,:B,:V,:K,:Tl,:Yb,:Sb,:Sn,:Ag,:Pd,:Co,:Se,:Ti,:Zn,:H,:Li,:Ge,:Cu,:Au,:Ni,:Cd,:In,:Mn,:Zr,:Cr,:Pt,:Hg,:Pb,:Unknown],
    "degree" => [0,1,2,3,4,5,6,7,8,9,10],
    "implicit_valence" => [0,1,2,3,4,5,6],
    "formal_charge" => [],#[-1,0,1,2], # indicates using real value
    "radical_electrons" => [],#[0,1,2,3], # indicates using use real value
    "hybridization" => [chem.rdchem.HybridizationType.SP, chem.rdchem.HybridizationType.SP2,chem.rdchem.HybridizationType.SP3,chem.rdchem.HybridizationType.SP3D, chem.rdchem.HybridizationType.SP3D2],#[2,3,4,5,6],
    "aromaticity" => [],#[0,1], # indicates using real value
    "total_H_num" => [0,1,2,3,4]
)

# bins for bond features
const bond_feat_bins = Dict(
    "bond_type" =>[chem.rdchem.BondType.SINGLE, chem.rdchem.BondType.DOUBLE,chem.rdchem.BondType.TRIPLE, chem.rdchem.BondType.AROMATIC],#[1,2,3,12],
    "isConjugated" => [],# [0,1],# indicates using real value
    "isInring" => []#[0,1],# indicates using real value
)


"""
input:
    atom (RDKit.Atom)
    feature_list(Array{String}): list of atomic features of interest
output:
    atom_feat (Dict): Dictionary of atomic features
"""
function atom_features(atom; feature_list::Array{String}=atom_feat_list)
    atom_feat = Dict()
    for f in feature_list
        @assert f in atom_feat_list "Feature "*f*" is not in supported atomic features list: [symbol,degree,implicit_valence,formal_charge,radical_electrons,hybridization,aromaticity,total_H_num]"
        val = rdkit_atom_function[f](atom)
        bins = atom_feat_bins[f]
        if isempty(bins) # use real value instead of one-hot
            atom_feat[f] = Float32(val)
        else
            onehot_vec = zeros(Float32,size(bins,1))
            try  # check if the value can be found in the given bins
                onehot_vec = Float32.(val .== bins)
            catch
                @warn f*" value out of range: "*string(val)
                onehot_vec = Float32.(bins[end] .== bins)
            end
            atom_feat[f] = onehot_vec
        end
    end
    return atom_feat
end

"""
input:
    bond (RDKit.Bond)
    feature_list(Array{String}): list of bond features of interest
output:
    bond_feat (Dict): Dictionary of bond features
"""
function bond_features(bond; feature_list::Array{String}=bond_feat_list)
    bond_feat = Dict()
    for f in feature_list
        @assert f in bond_feat_list "Feature "*f*" is not in supported bond features list: [bond_type, isConjugated, isInring]"
        val = rdkit_bond_function[f](bond)
        bins = bond_feat_bins[f]
        if isempty(bins) # use real value instead of one-hot
            bond_feat[f] = Float32(val)
        else
            onehot_vec = zeros(Float32,size(bins,1))
            try  # check if the value can be found in the given bins
                onehot_vec = Float32.(val .== bins)
            catch
                @warn f*" value out of range: "*string(val)
                onehot_vec = Float32.(bins[end] .== bins)
            end
            bond_feat[f] = onehot_vec
        end
    end
    return bond_feat
end

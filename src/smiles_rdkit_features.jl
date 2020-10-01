## Basic functions
using DataStructures
using Flux: onehot
using PyCall
chem = pyimport_conda("rdkit.Chem", "rdkit", "conda-forge")

# currently supported atomic features and bond features
const atom_feat_list = ["symbol","degree","implicit_valence","formal_charge","radical_electrons","hybridization","aromaticity","total_H_num" ]
const bond_feat_list = ["bond_type","isConjugated","isInring"]

# Map atomic feature names to RDKit functions
const rdkit_atom_functions = Dict(
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
const rdkit_bond_functions = Dict(
    "bond_type"=>chem.Bond.GetBondType,
    "isConjugated"=>chem.Bond.GetIsConjugated,
    "isInring"=>chem.Bond.IsInRing
)

# bins for atomic features
const atom_feat_bins = Dict(
    "symbol" => ["C","N","O","S","F","Si","P","Cl","Br","Mg","Na","Ca","Fe","As","Al","I","B","V","K","Tl","Yb","Sb","Sn","Ag","Pd","Co","Se","Ti","Zn",'H',"Li","Ge","Cu","Au","Ni","Cd","In","Mn","Zr","Cr","Pt","Hg","Pb","Unknown"],
    "degree" => [0,1,2,3,4,5,6,7,8,9,10],
    "implicit_valence" => [0,1,2,3,4,5,6],
    "formal_charge" => [nothing],  # [-1,0,1,2], # indicates using real value
    "radical_electrons" => [nothing],  # [0,1,2,3], # indicates using use real value
    "hybridization" => [chem.rdchem.HybridizationType.SP, chem.rdchem.HybridizationType.SP2,chem.rdchem.HybridizationType.SP3,chem.rdchem.HybridizationType.SP3D, chem.rdchem.HybridizationType.SP3D2],  # [2,3,4,5,6],
    "aromaticity" => [nothing],  # [0,1], # indicates using real value
    "total_H_num" => [0,1,2,3,4]
)

# bins for bond features
const bond_feat_bins = Dict(
    "bond_type" => [chem.rdchem.BondType.SINGLE, chem.rdchem.BondType.DOUBLE,chem.rdchem.BondType.TRIPLE, chem.rdchem.BondType.AROMATIC],  #[1,2,3,12],
    "isConjugated" => [nothing],  # [0,1],# indicates using real value
    "isInring" => [nothing]  #[0,1],# indicates using real value
)


"""
input:
    atom (RDKit.Atom)
    feature_list(Array{String}): list of atomic features of interest
output:
    atom_feat (OrderedDict): Dictionary of atomic features
"""
function smiles_atom_features(atom; feature_list::Array{String}=atom_feat_list)
    atom_feat = OrderedDict()
    for f in feature_list
        @assert f in atom_feat_list "Feature "*f*" is not in supported atomic features list: [symbol,degree,implicit_valence,formal_charge,radical_electrons,hybridization,aromaticity,total_H_num]"
        val = rdkit_atom_functions[f](atom)
        bins = atom_feat_bins[f]
        if isnothing(bins[1]) # use real value instead of one-hot
            atom_feat[f] = Float32(val)
        else
            onehot_vec = zeros(Float32,size(bins,1))
            try  # check if the value can be found in the given bins
                onehot_vec = Float32.(val .== bins)
            catch
                @warn f*" value out of range: "*string(val)*", use the max accpeted value instead"
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
    bond_feat (OrderedDict): Dictionary of bond features
"""
function smiles_bond_features(bond; feature_list::Array{String}=bond_feat_list)
    bond_feat = OrderedDict()
    for f in feature_list
        @assert f in bond_feat_list "Feature "*f*" is not in supported bond features list: [bond_type, isConjugated, isInring]"
        val = rdkit_bond_functions[f](bond)
        bins = bond_feat_bins[f]
        if isnothing(bins[1]) # use real value instead of one-hot
            bond_feat[f] = Float32(val)
        else
            onehot_vec = zeros(Float32,size(bins,1))
            try  # check if the value can be found in the given bins
                onehot_vec = Float32.(val .== bins)
            catch
                @warn f*" value out of range: "*string(val)*", use the max accpeted value instead"
                onehot_vec = Float32.(bins[end] .== bins)
            end
            bond_feat[f] = onehot_vec
        end
    end
    return bond_feat
end

"""
Pair features that will be used in weave_featurize()
input:
    mol (RDKit SMILES mol): mol coverted from SIMLES
    bond_feat(OrderedDict): list of bond features of mol
    adj_list(Array{Array{Int}}): adjacency list
    bf_len(Int): bond feature length
    max_distance(Int): cut-off distance
output:
    features (Array{Float32,2}): pair features for mol
"""
function pair_features(mol, bond_feats::OrderedDict, adj_list; bf_len = 6, max_distance = 7)
    N = mol.GetNumAtoms()  # number of atoms
    features = zeros(Float32,N,N,bf_len+max_distance+1) # declear feature matrix
    if length(bond_feats) == 0
        @warn "No bonds found in "*string(chem.MolToSmiles(mol))*"! Pair feature matrix filled with zeros!"
        return reshape(features,(N*N,bf_len+max_distance+1))
    end
    rings = mol.GetRingInfo().AtomRings()
    for a1 in 1:N
        for a2 in adj_list[a1]
            # first bf_len features are bond features
            index = a1<a2 ? (a1,a2) : (a2,a1)
            features[a1, a2, 1:bf_len] = bond_feats[index]
        end
        for ring in rings
            ring = ring .+ 1  # Julia index starts from 1
            if a1 in ring
                # `bf_len+1`-th feature is if the pair of atoms are in the same ring
                features[a1, collect(ring), bf_len+1] .= 1
                features[a1, a1, bf_len+1] = 0
            end
        end
        distance_matrix = zeros(Float32,N,max_distance)
        distances = _get_distance!(distance_matrix, a1, adj_list, max_distance)
        features[a1,:,bf_len+2:end] = distance_matrix
    end
    features = vcat([features[i,:,:] for i=1:size(features,1)]...)
    return features
end

"""
Helper function for pair_features()
Finds distance between the input atom and all other atoms using BFS

input:
    distance_matrix(Array{Int,2}): predefined matrix to be filled in with distances
    a1(Int): index for current atom
    adj_list(Array{Int}): adjacency list for a1
    max_distance(Int): cut-off distance
"""
function _get_distance!(distance_matrix,a1, adj_list, max_distance = 7)
    visited = falses(length(adj_list))
    visited[a1] = true
    cur_list = adj_list[a1]
    cur_distance = 1
    while length(cur_list)!=0 && cur_distance<=max_distance
        next_list = []
        for a2 in cur_list
            if visited[a2]
                continue
            end
            visited[a2] = true
            distance_matrix[a2,cur_distance] = 1
            for a3 in adj_list[a2]
                if ! visited[a3]
                    push!(next_list,a3)
                end
            end
        end
        cur_list = next_list
        cur_distance += 1
    end
end

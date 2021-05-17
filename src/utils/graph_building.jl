module GraphBuilding

export build_graph
export weights_voronoi, weights_cutoff
export inverse_square, exp_decay

using PyCall
using ChemistryFeaturization
using Serialization

# options for decay of bond weights with distance...
# user can of course write their own as well
inverse_square(x) = x^-2.0
exp_decay(x) = exp(-x)

# TODO: figure out best/idiomatic way to pass through the keyword arguments, surely the copy/paste is not it
# TODO: option to featurize here?
"""
Function to build graph from a file storing a crystal structure (currently supports anything ase.io.read can read in). Returns an AtomGraph object.

# Arguments
## Required Arguments
- `file_path::String`: Path to ASE-readable file containing a molecule/crystal structure

## Keyword Arguments
- `normalize_weights::Bool=true`: Whether to rescale graph weights such that the maximum value is 1.0 (recommended)
- `use_voronoi::bool`: if true, use Voronoi method for neighbor lists, if false use cutoff method

    (The rest of these parameters are only used if `use_voronoi==false`)

- `cutoff_radius::Real=8.0`: cutoff radius for atoms to be considered neighbors (in angstroms)
- `max_num_nbr::Integer=12`: maximum number of neighbors to include (even if more fall within cutoff radius)
- `dist_decay_func::Function=inverse_square`: function to determine falloff of graph edge weights with neighbor distance

"""
function build_graph(file_path::String; use_voronoi::Bool=false, cutoff_radius::Real=8.0, max_num_nbr::Integer=12, dist_decay_func::Function=inverse_square, normalize_weights::Bool=true)
    aseio = pyimport_conda("ase.io", "ase", "conda-forge")
    atoms_object = aseio.read(file_path)

    # list of atom symbols
    atom_ids = [get(atoms_object, i-1).symbol for i=1:length(atoms_object)]

    # check if any nonperiodic BC's
    nonpbc = any(.!atoms_object.pbc)
    local cant_voronoi = false
    if nonpbc & use_voronoi
        @warn "Voronoi edge weights are not supported if any direction in the structure is nonperiodic. Using cutoff weights method..."
        cant_voronoi = true
    end

    if use_voronoi && !cant_voronoi
        s = pyimport_conda("pymatgen.core.structure", "pymatgen", "conda-forge")
        pmgase = pyimport_conda("pymatgen.io.ase", "pymatgen", "conda-forge")
        aa = pmgase.AseAtomsAdaptor()
        struc = aa.get_structure(atoms_object)
        weight_mat = weights_voronoi(struc)
    else
        nl = pyimport_conda("ase.neighborlist", "ase", "conda-forge")
        is, js, dists = nl.neighbor_list("ijd", atoms_object, cutoff_radius)
        weight_mat = weights_cutoff(is.+1, js.+1, dists; max_num_nbr=max_num_nbr, dist_decay_func=dist_decay_func)
    end

    if normalize_weights
        weight_mat = weight_mat ./ maximum(weight_mat)
    end

    return weight_mat, atom_ids
end

"""
Build graph using neighbor number cutoff method adapted from original CGCNN. Note that `max_num_nbr` is a "soft" max, in that if there are more of the same distance as the last, all of those will be added.

# TODO
- option to cut off by nearest, next-nearest, etc. by DISTANCE rather than NUMBER of neighbors
"""
function weights_cutoff(is, js, dists; max_num_nbr=12, dist_decay_func=inverse_square)
    # sort by distance
    ijd = sort([t for t in zip(is, js, dists)], by = t->t[3])

    # initialize neighbor counts
    num_atoms = maximum(is)
    local nb_counts = Dict(i=>0 for i=1:num_atoms)
    local longest_dists = Dict(i=>0.0 for i in 1:num_atoms)

    # iterate over list of tuples to build edge weights...
    # note that neighbor list double counts so we only have to increment one counter per pair
    weight_mat = zeros(Float32, num_atoms, num_atoms)
    for (i,j,d) in ijd
        # if we're under the max OR if it's at the same distance as the previous one
        if nb_counts[i] < max_num_nbr || isapprox(longest_dists[i], d)
            weight_mat[i, j] += dist_decay_func(d)
            longest_dists[i] = d
            nb_counts[i] += 1
        end
    end

    # average across diagonal, just in case
    weight_mat = 0.5 .* (weight_mat .+ weight_mat')
end

"""
Build graph using neighbors from faces of Voronoi polyedra and weights from areas. Based on the approach from https://github.com/ulissigroup/uncertainty_benchmarking/blob/aabb407807e35b5fd6ad06b14b440609ae09e6ef/BNN/data_pyro.py#L268
"""
function weights_voronoi(struc)
    num_atoms = size(struc)[1]
    sa = pyimport_conda("pymatgen.analysis.structure_analyzer", "pymatgen", "conda-forge")
    vc = sa.VoronoiConnectivity(struc)
    conn = vc.connectivity_array
    weight_mat = zeros(Float32, num_atoms, num_atoms)
    # loop over central atoms
    for atom_ind in 1:size(conn)[1]
        # loop over neighbor atoms
        for nb_ind in 1:size(conn)[2]
            # loop over each possible PBC image for chosen image
            for image_ind in 1:size(conn)[3]
                # only add as neighbor if atom is not current center one AND there is connectivity to image
                if (atom_ind != image_ind) & (conn[atom_ind, nb_ind, image_ind] != 0)
                    weight_mat[atom_ind, nb_ind] += conn[atom_ind, nb_ind, image_ind]/maximum(conn[atom_ind, :, :])
                end
            end
        end
    end

    # average across diagonal (because neighborness isn't strictly symmetric in the way we're defining it here)
    weight_mat = 0.5.* (weight_mat .+ weight_mat')
end

# TODO: graphs from SMILES via OpenSMILES.jl

end
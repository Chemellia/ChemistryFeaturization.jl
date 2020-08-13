using PyCall

"""
TODO:
 [x] tabulate element features in a dataframe (just in Python for now)
 [ ] tabulate every oxidation state and associated ionic radii (also w/Python)
 [ ] PyCall using oxi_state_guesses from pymatgen inside featurization to figure these out, then pull radius from table
 [ ] maybe in future: Shannon radii (need a robust way to determine coordination)
"""

# options for decay of bond weights with distance...
inverse_square(x) = x^-2.0
exp_decay(x) = exp(-x)

"""
Function to build graph from a CIF file of a crystal structure. Returns a tuple of the adjacency matrix and the list of elements in order of the graph nodes.


# Arguments
- `cif_path::String`: path to CIF file
- `use_voronoi::bool`: if true, use Voronoi method for neighbor lists, if false use cutoff method

    (The rest of these parameters are only used if use_voronoi is false)

- `radius::Float=8.0`: cutoff radius for atoms to be considered neighbors (in angstroms)
- `max_num_nbr::Integer=12`: maximum number of neighbors to include (even if more fall within cutoff radius)
- `dist_decay_func`: function (e.g. inverse_square or exp_decay) to determine falloff of graph edge weights with neighbor distance
"""
function build_graph(cif_path; use_voronoi=true, radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square, normalize=true)
    s = pyimport("pymatgen.core.structure")
    c = s.Structure.from_file(cif_path)
    num_atoms = size(c)[1]

    # list of atom symbols
    atom_ids = [site_element(s) for s in c]

    if use_voronoi
        weight_mat = build_graph_voronoi(c, num_atoms)
    else
        weight_mat = build_graph_cutoff(c, radius, max_num_nbr, dist_decay_func)
    end

    # normalize weights
    if normalize
        weight_mat = weight_mat ./ maximum(weight_mat)
    end

    return (weight_mat, atom_ids)
end

"""
Build graph using neighbors from faces of Voronoi polyedra and weights from areas. Based on the approach from https://github.com/ulissigroup/uncertainty_benchmarking/blob/aabb407807e35b5fd6ad06b14b440609ae09e6ef/BNN/data_pyro.py#L268
"""
function build_graph_voronoi(struc)
    num_atoms = size(struc)[1]
    sa = pyimport("pymatgen.analysis.structure_analyzer")
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

"""
Build graph using neighbor number cutoff method adapted from original CGCNN. Note that `max_num_nbr` is a "soft" max, in that if there are more of the same distance as the last, all of those will be added.

# TODO
- option to cut off by nearest, next-nearest, etc. by DISTANCE rather than NUMBER of neighbors
"""
function build_graph_cutoff(struc; radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square)
    #= find neighbors, requires a cutoff radius
    returns a NxM Array of PyObject PeriodicSite
    ...except when it returns a list of N of lists of length M...
    each PeriodicSite is (site, distance, index, image)
    N is num sites in crystal, M is num neighbors
    =#
    all_nbrs = c.get_all_neighbors(radius, include_index=true)

    # sort by distance
    # returns list of length N of lists of length M
    if length(size(all_nbrs)) == 2
        all_nbrs = [sort(all_nbrs[i,:], lt=(x,y)->isless(site_distance(x), site_distance(y))) for i in 1:num_atoms]
    elseif length(size(all_nbrs)) == 1
        all_nbrs = [sort(all_nbrs[i][:], lt=(x,y)->isless(site_distance(x), site_distance(y))) for i in 1:num_atoms]
    end

    # iterate through each list of neighbors (corresponding to neighbors of a given atom) to find bonds (eventually, graph edges)
    weight_mat = zeros(Float32, num_atoms, num_atoms)
    for atom_ind in 1:num_atoms
        this_atom = get(c, atom_ind-1)
        atom_nbs = all_nbrs[atom_ind]
        # iterate over each neighbor...
        for nb_num in 1:size(all_nbrs[atom_ind])[1]
            nb = atom_nbs[nb_num]
            global nb_ind = site_index(nb)
            # if we're under the max, add it for sure
            if nb_num < max_num_nbr
                weight_mat[atom_ind, nb_ind] = weight_mat[atom_ind, nb_ind] + dist_decay_func(site_distance(nb))
            # if we're at/above the max, add if distance is the same
            else
                # check we're not on the last one
                if nb_ind < size(atom_nbs)[1] - 1
                    next_nb = atom_nbs[nb_ind + 1]
                    # add another bond if it's the exact same distance to the next neighbor in the list
                    if are_equidistant(nb, next_nb)
                        weight_mat[atom_ind, nb_ind] = weight_mat[atom_ind, nb_ind] + dist_decay_func(site_distance(nb))
                    end
                end
            end
        end
    end

    # average across diagonal (because neighborness isn't strictly symmetric in the way we're defining it here)
    weight_mat = 0.5.* (weight_mat .+ weight_mat')
end

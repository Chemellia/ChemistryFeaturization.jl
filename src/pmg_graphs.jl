#=
Building graphs from CIF files using PyCall to the pymatgen package.
=#

using PyCall
using ChemistryFeaturization
using Glob
using Serialization

# options for decay of bond weights with distance...
inverse_square(x) = x^-2.0
exp_decay(x) = exp(-x)

# a few helper fcns for readability...
"Return the index of a given site in the structure."
site_index(site) = convert(UInt16, get(site, 2)) + 1

"Return the distance associated with a site in a neighbor list."
site_distance(site) = convert(Float64, get(site, 1))

# these next two functions return the information for the first species in a list – there should only be one because otherwise the structure would be disordered and we probably shouldn't be building a graph...(or maybe we add some fancy functionality later to do superpositions of species?)

"Return atomic number associated with a site."
site_atno(site) = [e.Z for e in site.species.elements][1]

"Return atomic symbol associated with a site."
site_element(site) = [e.symbol for e in site.species.elements][1]


"""
    are_equidistant(site1, site2)

Check if site1 and site2 are equidistant to within tolerance atol, in angstroms (for cutting off neighbor lists consistently).

Note that this only works if site1 and site2 are from a neighbor list from the same central atom.
"""
function are_equidistant(site1, site2, atol=1e-4)
    isapprox(site_distance(site1), site_distance(site2), atol=atol)
end

# TODO: figure out best/idiomatic way to pass through the keyword arguments, surely the copy/paste is not it
"""
Function to build graph from a file storing a crystal structure (currently supports .cif files and ASE .traj files). Returns an AtomGraph object.

# Arguments
- `cif_path::String`: path to CIF file
- `use_voronoi::bool`: if true, use Voronoi method for neighbor lists, if false use cutoff method

    (The rest of these parameters are only used if use_voronoi is false)

- `radius::Float=8.0`: cutoff radius for atoms to be considered neighbors (in angstroms)
- `max_num_nbr::Integer=12`: maximum number of neighbors to include (even if more fall within cutoff radius)
- `dist_decay_func`: function (e.g. inverse_square or exp_decay) to determine falloff of graph edge weights with neighbor distance
"""
# TODO: featurize here
function build_graph(file_path; use_voronoi=false, radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square, normalize=true)
    s = pyimport("pymatgen.core.structure")

    # TODO: this bit can probably be abstracted out to another fcn...
    if file_path[end-3:end]==".cif"
        c = s.Structure.from_file(cif_path)
    elseif file_path[end-4:end]==".traj"
        aseio = pyimport("ase.io")
        pmgase = pyimport("pymatgen.io.ase")
        atoms_object = aseio.read(file_path)
        aa = pmgase.AseAtomsAdaptor()
        c = aa.get_structure(atoms_object)
    else
        error("file format not supported")
    end

    num_atoms = size(c)[1]

    # list of atom symbols
    atom_ids = [site_element(s) for s in c]

    if use_voronoi
        weight_mat = weights_voronoi(c)
    else
        weight_mat = weights_cutoff(c; radius=radius, max_num_nbr=max_num_nbr, dist_decay_func=dist_decay_func)
    end

    # normalize weights
    if normalize
        weight_mat = weight_mat ./ maximum(weight_mat)
    end

    g = SimpleWeightedGraph{Int32}(Float32.(weight_mat))
    return AtomGraph(g, atom_ids)
end

"""
Build graph using neighbors from faces of Voronoi polyedra and weights from areas. Based on the approach from https://github.com/ulissigroup/uncertainty_benchmarking/blob/aabb407807e35b5fd6ad06b14b440609ae09e6ef/BNN/data_pyro.py#L268
"""
function weights_voronoi(struc)
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
function weights_cutoff(struc; radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square)
    #= find neighbors, requires a cutoff radius
    returns a NxM Array of PyObject PeriodicSite
    ...except when it returns a list of N of lists of length M...
    each PeriodicSite is (site, distance, index, image)
    N is num sites in crystal, M is num neighbors
    =#
    num_atoms = size(struc)[1]
    all_nbrs = struc.get_all_neighbors(radius, include_index=true)

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
        this_atom = get(struc, atom_ind-1)
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


"""
Function to build and serialize to file a batch of CIFs, optionally featurizing them as well. Saved .jls files will have the same names as the .cif ones but with the extensions modified.

# Arguments
- `cif_folder::String`: path to folder containing CIF files
- `output_folder::String`: path to folder where .jls files containing AtomGraph objects should be saved (will be created if it doesn't exist already)

If you want to featurize the graphs, at least the vector of `AtomFeat` objects describing the featurization procedure is required, and optionally the Dict mapping from elemental symbols to feature vectors (if it is not provided, it will be generated).

Other optional arguments are the optional arguments to `build_graph`: `use_voronoi`, `radius`, `max_num_nbr`, `dist_decay_func`, `normalize`

This function does not return anything.
    TODO: decide if there should be an option to return the graphs
"""
function build_graphs_from_cifs(cif_folder::String, output_folder::String; atom_featurevecs=Dict{String, Vector{Float32}}(), featurization=AtomFeat[], use_voronoi=false, radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square, normalize=true)
    # check if input folder exists and contains CIFs, if not throw error
    ciflist = glob(joinpath(cif_folder, "*.cif"))
    if length(ciflist)==0
        error("No CIF's in provided CIF directory!")
    end

    # check if output folder exists, if not create it
    if !isdir(output_folder)
        mkdir(output_folder)
        @info "Output path provided did not exist, creating folder there."
    end

    # check if there is enough information to actually featurize
    featurize = length(featurization)>0
    atom_featurevecs = (length(atom_featurevecs)>0 & featurize) ? atom_featurevecs : make_feature_vectors(featurization)[1]

    if (length(atom_featurevecs)>0) & (length(featurization)==0)
        @warn "You have supplied only feature vectors but no featurization scheme, so graphs will be built but not featurized."
    end

    # loop over CIFs
    for cif in ciflist
        id = split(split(cif, "/")[end], ".")[1]
        ag = build_graph(cif; use_voronoi=use_voronoi, radius=radius, max_num_nbr=max_num_nbr, dist_decay_func=dist_decay_func, normalize=normalize)
        if featurize
            add_features!(ag, atom_featurevecs, featurization)
        end
        graph_path = joinpath(output_folder, string(id, ".jls"))
        serialize(graph_path, ag)
    end
end

# alternate call signature where featurization is generated
function build_graphs_from_cifs(cif_folder::String, output_folder::String; feature_names::Vector{Symbol}, nbins::Vector{<:Integer}=default_nbins*ones(Int64, size(feature_names,1)), logspaced=false)
    atom_featurevecs, featurization = make_feature_vectors(build_atom_feats(feature_names; nbins=nbins, logspaced=logspaced))
    build_graphs_from_cifs(cif_folder, output_folder; featurization=featurization, atom_featurevecs = atom_featurevecs)
end


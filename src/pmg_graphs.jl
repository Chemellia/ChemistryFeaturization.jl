#=
Building graphs from CIF files using PyCall to the pymatgen package.
=#
module graph_building
export weights_voronoi, weights_cutoff

using SimpleWeightedGraphs
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
- `struc`: pymatgen Structure object
- `use_voronoi::bool`: if true, use Voronoi method for neighbor lists, if false use cutoff method

    (The rest of these parameters are only used if use_voronoi is false)

- `radius::Float=8.0`: cutoff radius for atoms to be considered neighbors (in angstroms)
- `max_num_nbr::Integer=12`: maximum number of neighbors to include (even if more fall within cutoff radius)
- `dist_decay_func`: function (e.g. inverse_square or exp_decay) to determine falloff of graph edge weights with neighbor distance
"""
function build_graph(file_path::String; use_voronoi=false, radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square, normalize=true)
    # see if this fixes issues on Windows
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
        is, js, dists = nl.neighbor_list("ijd", atoms_object, radius)
        weight_mat = weights_cutoff(is.+1, js.+1, dists; max_num_nbr=max_num_nbr, dist_decay_func=dist_decay_func)
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

# TODO: might want a version that only serializes for cases where full array of AtomGraphs won't fit in memory?
"""
Function to build and, optionally, serialize to file a batch of input files, optionally featurizing them as well. Saved .jls files will have the same names as the .cif (or .traj, or .xyz, or ...) ones but with the extensions modified.

# Arguments
- `input_folder::String`: path to folder containing CIF files
- `output_folder::String` (optional): path to folder where .jls files containing AtomGraph objects should be saved (will be created if it doesn't exist already)

If you want to featurize the graphs, at least the vector of `AtomFeat` objects describing the featurization procedure is required, and optionally the Dict mapping from elemental symbols to feature vectors (if it is not provided, it will be generated).

Other optional arguments are the optional arguments to `build_graph`: `use_voronoi`, `radius`, `max_num_nbr`, `dist_decay_func`, `normalize`
"""
function build_graphs_batch(input_folder::String, featurization=AtomFeat[]; atom_featurevecs=Dict{String, Vector{Float32}}(), use_voronoi=false, radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square, normalize=true, output_folder="")
    # check if input folder exists and contains things, if not throw error
    file_list = readdir(input_folder, join=true)
    length(file_list)!=0 || throw(ArgumentError("No files in input directory!"))

    # check if output folder exists, if not create it
    to_serialize = false
    if output_folder != ""
        to_serialize=true
        if !isdir(output_folder)
            mkpath(output_folder)
            @info "Output path provided did not exist, creating folder there."
        end
    end

    # check if there is enough information to actually featurize
    to_featurize = length(featurization)>0
    atom_featurevecs = (length(atom_featurevecs)>0 & to_featurize) ? atom_featurevecs : make_feature_vectors(featurization)[1]

    if (length(atom_featurevecs)>0) & (length(featurization)==0)
        @warn "You have supplied only feature vectors but no featurization scheme, so graphs will be built but not featurized."
    end

    # loop over CIFs and build graphs from all files that we can...
    graphs = AtomGraph[]
    for file in file_list
        id = split(splitpath(file)[end], ".")[1]
        local ag
        try 
            ag = build_graph(file; use_voronoi=use_voronoi, radius=radius, max_num_nbr=max_num_nbr, dist_decay_func=dist_decay_func, normalize=normalize)
        catch
            @warn "Unable to build graph for $file"
            continue
        end
        if to_featurize
            add_features!(ag, atom_featurevecs, featurization)
        end
        if to_serialize
            graph_path = joinpath(output_folder, string(id, ".jls"))
            serialize(graph_path, ag)
        end
        push!(graphs, ag)
    end
    return graphs
end

# alternate call signature where featurization is generated
function build_graphs_batch(cif_folder::String, feature_names::Vector{Symbol}; nbins::Vector{<:Integer}=default_nbins*ones(Int64, size(feature_names,1)), logspaced=false, output_folder="")
    atom_featurevecs, featurization = make_feature_vectors(build_atom_feats(feature_names; nbins=nbins, logspaced=logspaced))
    build_graphs_batch(cif_folder, output_folder; featurization=featurization, atom_featurevecs = atom_featurevecs)
end

# function to read in graphs generated by `build_graphs_batch` to an Array of AtomGraph objects
function read_graphs_batch(graph_folder::String)
    files = [fp for fp in readdir(graph_folder, join=true) if fp[end-3:end]==".jls"]
    local graphs = AtomGraph[]
    for fpath in files
        push!(graphs, deserialize(fpath))
    end
    return graphs
end

end

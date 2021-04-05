#=
NB this code won't actually run because several variables don't ever get defined, it's just to show what usage looks like
=#

features = ["Group", "Row", "Block", "Atomic mass", "Atomic radius", "X"]
num_bins = [18, 9, 4, 16, 10, 10]
logspaced = [false, false, false, true, true, false]

atom_feature_vecs, featurization = make_feature_vectors(Symbol.(features), nbins=num_bins, logspaced=logspaced)

inputs = AtomGraph[]
# loop to deserialize precomputed AtomGraphs...
for r in eachrow(info)
    graph_path = joinpath(graph_dir, string("struc", lpad(r.Column1, 6, "0"), ".jls"))
    gr = deserialize(graph_path)
    feature_mat = hcat([atom_feature_vecs[e] for e in gr.elements]...)
    add_features!(gr, feature_mat, featurization)
    push!(inputs, gr)
end
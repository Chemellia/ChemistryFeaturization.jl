# TODO: proper docstring
#=
Featurization for `AtomGraph` objects that featurizes graph nodes only.
=#
struct GraphNodeFeaturization <: AbstractFeaturization
    atom_features::Vector{AtomFeature}
end

# TODO: function to compute total vector length

# docstring
function GraphNodeFeaturization(
    feature_names::Vector{String};
    nbins = default_nbins,
    logspaced = getindex.(Ref(default_log), feature_names),
)
    afs = map(zip(feature_names, nbins, logspaced)) do args
        AtomFeature(args[1], nbins = args[2], logspaced = args[3])
    end
    GraphNodeFeaturization(afs)
end
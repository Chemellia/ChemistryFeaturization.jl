#= FEATURIZATION OBJECTS
All such objects should define at least one list of <:AbstractFeature objects and either work according to the generic featurize! defined herein or dispatch featurize! if customized behavior is needed.
=#

module Featurizations

using Functors

# include...
include("graphnodefeaturization.jl")
include("weavefeaturization.jl")

# export...



abstract type AbstractFeaturization end

# generic featurize!
# this assumes that `a` has fields with names corresponding to each field in `fzn`, if not you need to dispatch this function to your specific case
# TODO: maybe add option to exclude field names from iteration over fzn?
# docstring
function featurize!(a<:AbstractAtoms, fzn<:AbstractFeaturization)
    # loop over fields in featurization, each one is a list of features
    # encode each feature in that list and assign the results to the
    # field of the same name in `a`
    for feats_list in fieldnames(fzn)
        encoded = reduce(vcat, map((x)->x(a), feats_list))
        setproperty!(a, feats_list, encoded)
    end
    a.featurization = fzn
end

end
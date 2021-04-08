module FeatureObjects

# link to guidance in docs about how to implement new feature types

# include...
include("atomfeat.jl")
include("pairfeat.jl")

# export...


abstract type AbstractFeature{Tn,Te} end

# magical encoding
function (f<:AbstractFeature{Tn,Te})(a<:AbstractAtoms)
    f.encode_f(a)
end

end

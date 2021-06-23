module Atoms

import ..ChemistryFeaturization.decode
import ..ChemistryFeaturization.AbstractType.AbstractAtoms

decode(a::AbstractAtoms) = decode(a.featurization, a.encoded_features)

export decode

include("atomgraph.jl")
export AtomGraph, visualize

end

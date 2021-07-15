module Atoms

using ..ChemistryFeaturization.AbstractType: AbstractAtoms

import ..ChemistryFeaturization.elements
elements(a::AbstractAtoms) = throw(MethodError(elements, a))
export elements

include("atomgraph.jl")
export AtomGraph, visualize, elements

end

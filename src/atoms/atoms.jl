"""
    Atoms

Atoms objects characterize structural representations of molecular structures.
These representations can either be standardized forms
(like graph representation using an adjacency matrix, point cloud representation, etc),
or be custom user-defined representations too.
"""
module Atoms

include("atomgraph.jl")
export AtomGraph, visualize

end

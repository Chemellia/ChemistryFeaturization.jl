module SpeciesFeatureUtils

using MolecularGraph

# some convenience SFD constructors mapping names of species features to MolecularGraph functions...

mg_elements = [d["Symbol"] for d in MolecularGraph.ATOMTABLE]

export sfd_names_props
sfd_names_props = Dict(
    "isaromatic" => Dict(
        :A => GraphMol,
        :compute_f => isaromatic,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [true, false],
    ),
)

end

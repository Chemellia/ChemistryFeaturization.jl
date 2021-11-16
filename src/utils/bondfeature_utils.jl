module BondFeatureUtils

using MolecularGraph

# NB this is duplicated from SpeciesFeatureUtils and thus should probably eventually be moved somewhere else
const mg_elements = [d["Symbol"] for d in MolecularGraph.ATOMTABLE]

export bfd_names_props
const bfd_names_props = Dict(
    "bondorder" => Dict(
        :A => GraphMol,
        :compute_f => bondorder,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [1, 2, 3],
    ),
    "isaromaticbond" => Dict(
        :A => GraphMol,
        :compute_f => isaromaticbond,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [true, false],
    ),
    "isringbond" => Dict(
        :A => GraphMol,
        :compute_f => isringbond,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [true, false],
    ),
    "isrotatable" => Dict(
        :A => GraphMol,
        :compute_f => isrotatable,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [true, false],
    ),
)

end
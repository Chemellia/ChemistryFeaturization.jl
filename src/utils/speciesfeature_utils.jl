module SpeciesFeatureUtils

using ...ChemistryFeaturization.Codec: DirectCodec

using MolecularGraph

# some convenience SFD constructors mapping names of species features to MolecularGraph functions...

const mg_elements = [d["Symbol"] for d in MolecularGraph.ATOMTABLE]

export sfd_names_props
const sfd_names_props = Dict(
    "hybridization" => Dict(
        :A => GraphMol,
        :compute_f => hybridization,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [:sp3, :sp2, :sp, :none],
    ),
    "hydrogenconnected" => Dict(
        :A => GraphMol,
        :compute_f => hydrogenconnected,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [0, 1, 2, 3], # I'm hoping this is right...
    ),
    "isaromatic" => Dict(
        :A => GraphMol,
        :compute_f => isaromatic,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [true, false],
    ),
    "isringatom" => Dict(
        :A => GraphMol,
        :compute_f => isringatom,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [true, false],
    ),
    "lonepair" => Dict(
        :A => GraphMol,
        :compute_f => lonepair,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [-1, 0, 1, 2, 3],
    ),
    "pielectron" => Dict(
        :A => GraphMol,
        :compute_f => pielectron,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [0, 1, 2], # also not certain this is correct
    ),
    "degree" => Dict(
        :A => GraphMol,
        :compute_f => nodedegree,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => collect(0:10),
    ),
    "radical_electrons" => Dict(
        :A => GraphMol,
        :compute_f => multiplicity,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [1, 2, 3],
    ),
    "multiplicity" => Dict(
        :A => GraphMol,
        :compute_f => multiplicity,
        :categorical => true,
        :encodable_elements => mg_elements,
        :possible_vals => [1, 2, 3],
    ),
    "implicithconnected" => Dict(
        :A => GraphMol,
        :compute_f => implicithconnected,
        :categorical => false,
        :encodable_elements => mg_elements,
        :codec => DirectCodec,
    ),
    "charge" => Dict(
        :A => GraphMol,
        :compute_f => charge,
        :categorical => false,
        :encodable_elements => mg_elements,
        :codec => DirectCodec,
    ),
)

end

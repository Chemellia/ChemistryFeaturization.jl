module ElementFeature

using DataFrames

include("elementfeature_utils.jl")
using ..ChemistryFeaturization.Data: atom_data_df
import ..ChemistryFeaturization: encodable_elements, get_value, default_codec, AbstractAtomFeatureDescriptor, OneHotOneCold

export ElementFeatureDescriptor, get_value, default_codec, encodable_elements

"""
    ElementFeatureDescriptor

A descriptor for features associated with individual atoms that depend only upon their elemental identity (and whose values can hence be determined from a lookup table).

## Fields
- `name::String`: Name of the feature
- `lookup_table::DataFrame`: table containing values of feature for every encodable element
"""
struct ElementFeatureDescriptor<: AbstractAtomFeatureDescriptor
    name::String
    lookup_table::DataFrame
    function ElementFeatureDescriptor(
        feature_name::String,
        lookup_table = atom_data_df,
    )
        colnames = names(lookup_table)
        @assert feature_name in colnames && "Symbol" in colnames "Your lookup table must have a column called :Symbol and one with the same name as your feature to be usable!"

        lookup_table = lookup_table[:, ["Symbol", feature_name]]
        dropmissing!(lookup_table)

        new(
            feature_name,
            lookup_table,
        )
    end
end

encodable_elements(efd::ElementFeatureDescriptor) = efd.lookup_table[:, :Symbol]

function get_value(efd::ElementFeatureDescriptor, a)
    @assert all([el in encodable_elements(efd) for el in elements(a)]) "Feature $(efd.name) cannot encode some element(s) in this structure!"

    feature_vals = efd.lookup_table[:, [:Symbol, Symbol(efd.name)]]
    map(
        el ->
            getproperty(feature_vals[feature_vals.Symbol.==el, :][1, :], Symbol(efd.name)),
        elements(a),
    )
end

default_codec(efd::ElementFeatureDescriptor) =
    OneHotOneCold(default_categorical(efd.name, efd.lookup_table), get_bins(efd.name, efd.lookup_table))

end
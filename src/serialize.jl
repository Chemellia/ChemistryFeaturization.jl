using ..ChemistryFeaturization.AbstractType: AbstractAtoms, AbstractFeaturization

"""
    SerializableEncodedFeature

Helps in serialization, and keeps a record of the `featurization` applied to an `atoms`, and the `encoded_features`
that are generated as a result.

Note: `encoded_features` will NOT change for a given atoms-featurization pair.
"""
struct SerializableEncodedFeature
    atoms::AbstractAtoms
    featurization::AbstractFeaturization
    encoded_features::Any

    # inner constructor - generates `encoded_features`
    SerializableEncodedFeature(atoms::AbstractAtoms, featurization::AbstractFeaturization) = new(atoms, featurization, featurize!(atoms, featurization))
end

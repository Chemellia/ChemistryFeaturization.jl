using ..ChemistryFeaturization.AbstractType: AbstractAtoms, AbstractFeaturization

"""
    FeaturizedAtoms

Container object for an `Atoms` object, a `featurization`, and the resulting
`encoded_features` from applying the `featurization` to the `atoms`.

Note: `encoded_features` will NOT change for a given atoms-featurization pair.
"""
struct FeaturizedAtoms
    atoms::AbstractAtoms
    featurization::AbstractFeaturization
    encoded_features::Any

    # inner constructor - generates `encoded_features`
    FeaturizedAtoms(atoms::AbstractAtoms, featurization::AbstractFeaturization) =
        new(atoms, featurization, featurize(atoms, featurization))
end

function decode(featurized_atoms::FeaturizedAtoms)
    decoded = decode(featurized_atoms.featurization, featurized_atoms.encoded_features)
    # is this loop universally applicable?
    # for (k, v) in decoded
    #     v["Symbol"] = sef.atoms.elements[k]
    # end
    return decoded
end

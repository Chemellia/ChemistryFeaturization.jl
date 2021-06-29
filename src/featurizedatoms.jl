using ..ChemistryFeaturization.AbstractType: AbstractAtoms, AbstractFeaturization

"""
    FeaturizedAtoms

Container object for an `Atoms` object, a `featurization`, and the resulting
`encoded_features` from applying the `featurization` to the `atoms`.

Note: `encoded_features` will NOT change for a given atoms-featurization pair.
"""
struct FeaturizedAtoms{A<:AbstractAtoms,F<:AbstractFeaturization}
    atoms::A
    featurization::F
    encoded_features::Any
    FeaturizedAtoms{A,F}(
        atoms,
        featurization,
    ) where {A<:AbstractAtoms,F<:AbstractFeaturization} =
        new(atoms, featurization, featurize(atoms, featurization))
end

FeaturizedAtoms(
    atoms::A,
    featurization::F,
) where {A<:AbstractAtoms,F<:AbstractFeaturization} =
    FeaturizedAtoms{A,F}(atoms, featurization)

function decode(featurized_atoms::FeaturizedAtoms)
    decoded = decode(featurized_atoms.featurization, featurized_atoms.encoded_features)
    # is this loop universally applicable?
    # for (k, v) in decoded
    #     v["Symbol"] = sef.atoms.elements[k]
    # end
    return decoded
end

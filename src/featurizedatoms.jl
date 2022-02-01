"""
    FeaturizedAtoms

Container object for an atomic structure object, a featurization, and the resulting
`encoded_features` from applying the `featurization` to the `atoms`.

## Fields
- `atoms`: object to be featurized
- `featurization`: Featurization scheme meant to be used for featurizing `atoms`
- `encoded_features`: The result of featurizing `atoms` using `featurization`

!!! note
    `encoded_features` will NOT change for a given atoms-featurization pair.
"""
struct FeaturizedAtoms{A,F<:AbstractFeaturization}
    atoms::A
    featurization::F
    encoded_features::Any
    FeaturizedAtoms{A,F}(atoms, featurization) where {A,F<:AbstractFeaturization} =
        new(atoms, featurization, encode(atoms, featurization))
end

FeaturizedAtoms(atoms::A, featurization::F) where {A,F<:AbstractFeaturization} =
    FeaturizedAtoms{A,F}(atoms, featurization)

# pretty printing, short
Base.show(io::IO, fa::FeaturizedAtoms) =
    print(io, string(typeof(fa), " with ", length(features(fa.featurization)), " features"))

# # pretty printing, long
function Base.show(io::IO, ::MIME"text/plain", fa::FeaturizedAtoms)
    st = string("FeaturizedAtoms:\n")
    st = string(st, "\tAtoms: $(fa.atoms)\n\tFeaturization: $(fa.featurization)")
    print(io, st)
end

"""
    decode(featurized_atoms::FeaturizedAtoms)

Decode a `FeaturizedAtoms` object, and return the decoded value.
"""
decode(featurized_atoms::FeaturizedAtoms) =
    decode(featurized_atoms.encoded_features, featurized_atoms.featurization)


"""
    featurize(atoms, featurization::AbstractFeaturization)

Featurize an `atoms` object using a `featurization` and return the
`FeaturizedAtoms` object created.
"""
featurize(atoms, featurization::AbstractFeaturization) =
    FeaturizedAtoms(atoms, featurization)

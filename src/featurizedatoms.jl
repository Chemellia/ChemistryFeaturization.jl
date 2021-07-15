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
        new(atoms, featurization, encode(featurization, atoms))
end

FeaturizedAtoms(
    atoms::A,
    featurization::F,
) where {A<:AbstractAtoms,F<:AbstractFeaturization} =
    FeaturizedAtoms{A,F}(atoms, featurization)

# pretty printing, short
function Base.show(io::IO, fa::FeaturizedAtoms)
    a, b = size(fa.encoded_features)
    st = string(typeof(fa), " with $a x $b encoded features")
    print(io, st)
end

# pretty printing, long
function Base.show(io::IO, ::MIME"text/plain", fa::FeaturizedAtoms)
    a, b = size(fa.encoded_features)
    st = string(typeof(fa), " with $a x $b encoded features:\n")
    st = string(st, "\tAtoms: $(fa.atoms)\n\tFeaturization: $(fa.featurization)")
    print(io, st)
end

"""
    decode(featurized_atoms::FeaturizedAtoms)

Decode a [FeaturizedAtoms](@ref) object, and return the decoded value.
"""
function decode(featurized_atoms::FeaturizedAtoms)
    decoded = decode(featurized_atoms.featurization, featurized_atoms.encoded_features)
    # is this loop universally applicable?
    # for (k, v) in decoded
    #     v["Symbol"] = sef.atoms.elements[k]
    # end
    return decoded
end


"""
    featurize(atoms::AbstractAtoms, featurization::AbstractFeaturization)

Featurize an `atoms` object using a `featurization` and return the
[FeaturizedAtoms](@ref) object created.
"""
featurize(atoms::AbstractAtoms, featurization::AbstractFeaturization) =
    FeaturizedAtoms(atoms, featurization)

export DummyFeaturization, featurize!
export encodable_elements, decode

using ..ChemistryFeaturization.AbstractType: AbstractFeaturization

"""
    DummyFeaturization

A dummy featurization that cannot actually encode or decode anything. For use primarily for populating output of model layers (e.g. AtomicGraphnets' AGNConv) that return an Atoms object, but with the encoded features transformed such that the original featurization is no longer applicable/valid.
"""
struct DummyFeaturization <: AbstractFeaturization end

encodable_elements(df::DummyFeaturization) = []

decode(fzn::DummyFeaturization, encoded_feature) = throw(ArgumentError("This featurization is just a dummy, likely created by a model layer such as `AGNConv`, and cannot actually decode encoded features."))

featurize!(a::AbstractAtoms, df::DummyFeaturization) = throw(ArgumentError("This featurization is just a dummy, likely created by a model layer such as `AGNConv`, and cannot actually encode features."))
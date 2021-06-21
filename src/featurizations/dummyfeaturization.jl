export DummyFeaturization, featurize!
export encodable_elements, decode

using ..ChemistryFeaturization.AbstractType: AbstractFeaturization

"""
    DummyFeaturization

A dummy featurization that cannot actually encode or decode anything. For use primarily for model layers (e.g. AtomicGraphnets' AGNConv) that return an Atoms object, but with the encoded features transformed such that the original featurization is no longer applicable/valid.
"""
struct DummyFeaturization <: AbstractFeaturization end

encodable_elements(df::DummyFeaturization) = []
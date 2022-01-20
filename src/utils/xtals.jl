using StaticArrays
using Xtals
using Zygote: @adjoint
using Zygote.ChainRulesCore
using NearestNeighbors
using LinearAlgebra

function rl(f_to_c::Array{Float64, 2})
   # the unit cell vectors are the columns of f_to_c
   a₁ = f_to_c[:, 1]
   a₂ = f_to_c[:, 2]
   a₃ = f_to_c[:, 3]

   # r = zeros(Float64, 3, 3)
   l1 = 2 * π * cross(a₂, a₃) / dot(a₁, cross(a₂, a₃))
   l2 = 2 * π * cross(a₃, a₁) / dot(a₂, cross(a₃, a₁))
   l3 = 2 * π * cross(a₁, a₂) / dot(a₃, cross(a₁, a₂))
   vcat(l1', l2', l3')
end

@adjoint function Xtals.reciprocal_lattice(x)
  Zygote.pullback(rl, x)
end

function replicate2(crystal::Crystal, repfactors::Tuple{Int, Int, Int})
    if Xtals.ne(crystal.bonds) != 0
        error("the crystal " * crystal.name * " has assigned bonds. to replicate, remove
        its bonds with `remove_bonds!(crystal)`. then use `infer_bonds(crystal)` to
        reassign the bonds")
    end

    assert_P1_symmetry(crystal)

    n_atoms = crystal.atoms.n * prod(repfactors)
    n_charges = crystal.charges.n * prod(repfactors)

    box = replicate(crystal.box, repfactors)

    xf_shift = Zygote.ignore() do
      rf = range.(0, repfactors .- 1, step = 1)
      x = repeat(collect.(sort(vec(collect(Iterators.product(rf...))))), inner = crystal.atoms.n)
      reduce(hcat, x)
    end

    # Repeat Atoms
    xf_raw = repeat(crystal.atoms.coords.xf, 1, prod(repfactors)) .+ xf_shift
    xf = xf_raw ./ repfactors
    
    frac = Frac(xf)
    species = repeat(crystal.atoms.species, inner = prod(repfactors))
    atoms = Xtals.Atoms(length(species), species, frac)

    # Repeat Charges
    q = repeat(crystal.charges.q, inner = prod(repfactors))
    xf_raw = repeat(crystal.charges.coords.xf, 1, prod(repfactors))
    xf = length(xf_raw) > 0 ? (xf_raw .+ xf_shift) ./ repfactors : xf_raw
    charges = Xtals.Charges(length(q), q, Frac(xf))
    
    return Crystal(crystal.name, box, atoms, charges, Xtals.MetaGraph(n_atoms), crystal.symmetry)
end

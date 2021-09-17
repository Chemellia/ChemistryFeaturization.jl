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

    charges = Xtals.Charges{Xtals.Frac}(n_charges)
    # atoms = Xtals.Atoms{Xtals.Frac}(n_atoms)

    xf_shift = Zygote.ignore() do
      x = repeat(collect.(sort(vec(collect(Iterators.product(0:2, 0:2, 0:2))))), inner = crystal.atoms.n)
      reduce(hcat, x)
    end

    xf_raw = repeat(crystal.atoms.coords.xf, 1, prod(repfactors)) .+ xf_shift
    xf = xf_raw ./ repfactors
    
    frac = Frac(xf)
    species = repeat(crystal.atoms.species, inner = prod(repfactors))
    atoms = Xtals.Atoms(length(species), species, frac)

    # charge_counter = 0
    # atom_counter = 0
    # for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
    #     xf_shift = 1.0 * [ra, rb, rc]

    #     # replicate atoms
    #     for i = 1:crystal.atoms.n
    #         atom_counter += 1

    #         atoms.species[atom_counter] = crystal.atoms.species[i]

    #         xf = crystal.atoms.coords.xf[:, i] + xf_shift
    #         atoms.coords.xf[:, atom_counter] = xf ./ repfactors
    #     end

    #     # replicate charges
    #     for i = 1:crystal.charges.n
    #         charge_counter += 1

    #         charges.q[charge_counter] = crystal.charges.q[i]

    #         xf = crystal.charges.coords.xf[:, i] + xf_shift
    #         charges.coords.xf[:, charge_counter] = xf ./ repfactors
    #     end
    # end
    # global reggatoms = atoms
    
    return Crystal(crystal.name, box, atoms, charges, Xtals.MetaGraph(n_atoms), crystal.symmetry)
end

function Xtals.distance(coords::Xtals.Frac, box::Xtals.Box, i, j, pbc)
  dxf = @views coords.xf[:, i] - coords.xf[:, j]
  norm.(eachcol(box.f_to_c * dxf))
end


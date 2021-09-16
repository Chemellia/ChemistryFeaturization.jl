using StaticArrays
using Xtals
using Zygote: @adjoint
using Zygote.ChainRulesCore

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

    xf_shift, ix = Zygote.ignore() do
      repeat(vec(collect.(Iterators.product(0:2, 0:2, 0:2))), inner = 27), repeat(1:crystal.atoms.n, outer = 27)
    end
    cols = map(ix) do i
      xf = @views crystal.atoms.coords.xf[:, i] + xf_shift[i]
      xf ./ repfactors
    end
    xf = reduce(hcat, cols)
    frac = Frac(xf)
    species = repeat(crystal.atoms.species, inner = 27)
    atoms = Xtals.Atoms(length(species), species, frac)

    # charge_counter = 0
    # for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
    #     xf_shift = 1.0 * [ra, rb, rc]

    #     # # replicate atoms
    #     # for i = 1:crystal.atoms.n
    #     #     atom_counter += 1

    #     #     atoms.species[atom_counter] = crystal.atoms.species[i]

    #     #     xf = crystal.atoms.coords.xf[:, i] + xf_shift
    #     #     atoms.coords.xf[:, atom_counter] = xf ./ repfactors
    #     # end

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

@adjoint function Xtals.distance_matrix(c::Xtals.Crystal, pbc)
  n = c.atoms.n
  Xtals.distance_matrix(c, pbc), Δ -> begin
  d = collect(Δ)
  dfc = similar(Δ, size(c.box.f_to_c))
  nt = Zygote.accum(Zygote.nt_nothing(c),
                    (atoms = Zygote.nt_nothing(c.atoms),
                     box = Zygote.nt_nothing(c.box)))
  for i = 1:n
    for j = (i+1):n
      y, back = Zygote.pullback(@view(c.atoms.coords.xf[:,i]),
                                @view(c.atoms.coords.xf[:,j]),
                                c.box.f_to_c) do x_, y_, fc
        # TODO: support for periodic boundary conditions
        # Xtals.distance(at, bo, i, j, pbc)
        dxf = x_ - y_
        norm(fc * dxf)
      end
      b = back(d[i,j])
      l = length(b[1])
      cix = vcat(repeat([CartesianIndex((i,i))], l),
                 repeat([CartesianIndex((i,j))], l))
      jix = vcat(repeat([CartesianIndex((i,i))], l),
                 repeat([CartesianIndex((j,i))], l))
      dfc = Zygote.accum(dfc, b[3])
      v = vcat(b[1], b[2])
      d[cix] = v
      d[jix] = v
    end
  end

  # Accumulation
  ant = (atoms = (coords = (xf = d,),),)
  bnt = (box = (f_to_c = dfc, ),)
  nt = Zygote.accum(nt, ant)
  nt = Zygote.accum(nt, bnt)
  (nt, nothing)
  end
end


function ChainRulesCore.rrule(::Type{SArray{D, T, ND, L}}, x...) where {D, T, ND, L}
  y = SArray{D, T, ND, L}(x...)
  function sarray_pb(Δy)
    Δy = map(t->eltype(x...)(t...), Δy)
    return NoTangent(), (Δy...,)
  end
  return y, sarray_pb
end


function index_works(crystal::Xtals.Crystal, n_atoms; cutoff_radius = 8.)
  tree = BruteTree(Cart(crystal.atoms.coords, crystal.box).x)

  is_raw = 13*n_atoms+1:14*n_atoms
  js_raw = inrange(tree,
                   Cart(crystal.atoms.coords[is_raw],
                        crystal.box).x,
                   cutoff_radius)

  split1 = map(zip(is_raw, js_raw)) do x
      [
          p for p in [(x[1], [j for j in js if j != x[1]]...) for js in x[2]] if
          length(p) == 2
      ]
  end
  ijraw_pairs = [(split1...)...]
end
Zygote.@nograd index_works

index_map(i, n_atoms) = (i - 1) % n_atoms + 1
Zygote.@nograd index_map

function neighbor_list2(crys::Crystal; cutoff_radius::Real = 8.0)
    n_atoms = crys.atoms.n

    # make 3 x 3 x 3 supercell and find indices of "middle" atoms
    # as well as index mapping from outer -> inner
    supercell = replicate2(crys, (3, 3, 3))

    # check for size of cutoff radius relative to size of cell
    min_celldim = min(crys.box.a, crys.box.b, crys.box.c)
    if cutoff_radius >= min_celldim
        @warn "Your cutoff radius is quite large relative to the size of your unit cell. This may cause issues with neighbor list generation, and will definitely cause a very dense graph. To avoid issues, I'm setting it to be approximately equal to the smallest unit cell dimension."
        cutoff_radius = 0.99 * min_celldim
    end


    ijraw_pairs = Zygote.ignore() do
      index_works(supercell, n_atoms, cutoff_radius = cutoff_radius)
    end
    dists = Xtals.distance_matrix(supercell, false)
    is = index_map.([t[1] for t in ijraw_pairs], n_atoms)
    js = index_map.([t[2] for t in ijraw_pairs], n_atoms)
    return is, js, dists
end

function build_graph2(
    crys::Crystal;
    cutoff_radius::Real = 8.0,
    max_num_nbr::Integer = 12,
    dist_decay_func::Function = inverse_square,
)

    is, js, dists = neighbor_list2(crys; cutoff_radius = cutoff_radius)
    weight_mat = weights_cutoff(
        is,
        js,
        dists;
        max_num_nbr = max_num_nbr,
        dist_decay_func = dist_decay_func,
    )
    return weight_mat, String.(crys.atoms.species)
end


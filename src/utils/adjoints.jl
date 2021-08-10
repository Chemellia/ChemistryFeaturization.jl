using Zygote # , ChainRulesCore
using Zygote: @adjoint
using LinearAlgebra

@adjoint function Base.Generator(f, iter)
  ys, backs = Zygote.unzip([Zygote.pullback(f, x) for x in iter])
  Base.Generator(f, iter), Δ -> begin
    b(d::Dict) = [back(v)[1] for (back,v) in zip(backs,values(d))]
    b(nt::NamedTuple{(:f, :iter)}) = [back(i)[1] for (back,i) in zip(backs,nt.iter)]
    (nothing, b(Δ))
  end
end

@adjoint function Base.Iterators.Zip(is)
  Zip_pullback(Δ) = (Zygote.unzip(Δ),)
  return Base.Iterators.Zip(is), Zip_pullback
end

@adjoint function Pair(a, b)
  Pair(a, b), Δ -> (Δ, nothing)
end

@adjoint function Dict(g::Base.Generator)
  ys, backs = Zygote.unzip([Zygote.pullback(g.f, args) for args in g.iter])
  Dict(ys...), Δ -> begin
    ∂d = first(backs)(Δ)[1]
    d = Dict(ys...)
    for (k,v) in pairs(d)
      d[k] = _zero(v)
    end
    (merge(d, ∂d), )
  end
end

_zero(x) = zero(x)
_zero(::Nothing) = nothing

@adjoint function _cutoff!(weight_mat, f, ijd,
                           nb_counts, longest_dists;
                           max_num_nbr = 12)
  y, ld = _cutoff!(weight_mat, f, ijd,
               nb_counts, longest_dists;
               max_num_nbr = max_num_nbr)
  function cutoff_pb((Δ,nt))
    Δ = collect(Δ)
    ld2 = Dict(keys(longest_dists) .=> zero.(values(longest_dists)))
    nc = Dict(keys(nb_counts) .=> zero.(values(nb_counts)))
    for (i,j,d) in ijd
      if nc[i] < max_num_nbr || isapprox(longest_dists[i], d)
        y_, back_ = Zygote.pullback(f, d)# weight_mat[i,j])
        Δ[i,j] *= back_(Δ[i,j])[1]
        ld2[i] += zero(d) # one(d)
        nc[i] += 0
      end
    end
    (Δ, nothing,
    collect(zip(fill(nothing, size(Δ,1)),
                fill(nothing, size(Δ,1)),
                diag(Δ))),
    nothing,
    nothing)
  end

  (y,ld), cutoff_pb
end
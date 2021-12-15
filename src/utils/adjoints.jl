using Zygote # , ChainRulesCore
using Zygote: @adjoint
using StaticArrays
using LinearAlgebra

@adjoint function Base.Iterators.Zip(is)
  Zip_pullback(Δ) = (Zygote.unzip(Δ.iter),)
  return Base.Iterators.Zip(is), Zip_pullback
end

@adjoint function Dict(g::Base.Generator)
  ys, backs = Zygote.unzip([Zygote.pullback(g.f, args) for args in g.iter])
  Dict(ys...), Δ -> begin
    dd = Dict(k => b(Δ)[1].second for (b,(k,v)) in zip(backs, pairs(Δ)))
    ((x for x in dd),)
  end
end

@adjoint function Base.Generator(f, args)
  Base.Generator(f, args), Δ -> (nothing, Δ)
end
@adjoint function Pair(k, v)
  Pair(k, v), Δ -> begin
    @show Δ
    (nothing, Δ[k])
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
    s = size(Δ)
    Δ = vec(collect(Δ))
    for (ix, (_,_,d)) in zip(eachindex(Δ), ijd)
      y_, back_ = Zygote.pullback(f, d)
      Δ[ix] *= back_(Δ[ix])[1]
    end
    (reshape(Δ, s), nothing,
    collect(zip(fill(nothing, size(Δ,1)),
                fill(nothing, size(Δ,1)),
                Δ)),
    nothing,
    nothing)
  end

  (y,ld), cutoff_pb
end

Zygote.@nograd Xtals.Charges{Xtals.Frac}

function Zygote.ChainRulesCore.rrule(::Type{SArray{D, T, ND, L}}, x...) where {D, T, ND, L}
  y = SArray{D, T, ND, L}(x...)
  function sarray_pb(Δy)
    Δy = map(t->eltype(x...)(t...), Δy)
    return NoTangent(), (Δy...,)
  end
  return y, sarray_pb
end

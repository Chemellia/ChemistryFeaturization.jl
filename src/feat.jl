struct Atom{T,A,B}
  f::T
  a::A
  b::B
  r::Bool
end

# Atom{f}(a,b) where f = Atom(a, b)

Atom(f, a, b) = Atom(f, a, b, true)
featurize(x, args...) = encode(x, args...)

# users can define this to scale with fields
# depending on features they are interested in
struct WeaveFeat
  atom_feats
  pair_feats
  bond_feats
  r::Bool
end

Functors.@functor WeaveFeat

function featurize(a::WeaveFeat, args...)
  Functors.fmap(a) do x
    # @show x
    featurize(x, args...)
  end
end

encode(a::Atom) = a.f(a) # Atom(:feat1, 2a.a, a.b)
encode(x) = x

feat1(a) = Atom(a.f, 2a.a, a.b) # we have the ones we are interested in
feat2(a) = Atom(a.f, a.a, 2a.b) # users can define more

feat1(a) = Flux.onehotbatch(a.a, 1:a.b)
feat2(a) = Flux.onehotbatch(a.a, 1:2a.b)

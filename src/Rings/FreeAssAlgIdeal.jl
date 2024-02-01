###############################################################################
#
#  Ideals of a free associative algebra
#
###############################################################################

# Currently ideal membership relies entirely on Singular, where a degree bound
# is imposed and a inconclusive answer may be returned. We can later add the
# Groebner machinery operating purely on the Oscar types, and hence not
# necessarily be confined to a degree bound.

@doc raw"""
    mutable struct FreeAssAlgIdeal{T} <: FreeAssAlgIdeal{T}

Two-sided ideal of a free associative algebra with elements of type `T`.
"""
mutable struct FreeAssAlgIdeal{T}
  gens::IdealGens{T}
  gb::IdealGens{T}

  function FreeAssAlgIdeal(R::FreeAssAlgebra, g::Vector{T}) where {T <: FreeAssAlgElem}
    r = new{T}()
    r.gens = IdealGens(R, g)
    return r
  end
  function FreeAssAlgIdeals(Ox::T, s::Singular.sideal) where {T <: FreeAssAlgElem}
    r = new{elem_type(T)}()
    r.gens = IdealGens(Ox, s)
    r.deg_bound = -1
    if s.isGB
      r.gb = r.gens
      r.deg_bound = base_ring(s).deg_bound
    end
    return r
  end
end

function AbstractAlgebra.expressify(a::FreeAssAlgIdeal; context = nothing)
  return Expr(:call, :ideal, [expressify(g, context = context) for g in collect(a.gens)]...)
end
@enable_all_show_via_expressify FreeAssAlgIdeal

@doc raw"""
    ideal(R::FreeAssAlgebra, g::Vector{T}) where T <: FreeAssAlgElem

Return the two-sided ideal of $R$ generated by $g$.
"""
function ideal(R::FreeAssAlgebra, g::Vector{T}) where T <: FreeAssAlgElem
  @assert all(x -> parent(x) == R, g) "parent mismatch"
  return FreeAssAlgIdeal(R, g)
end


function ideal(g::Vector{T}) where T <: FreeAssAlgElem
  @assert length(g) > 0 "cannot infer base ring"
  return MPolyIdeal(g)
end


function base_ring(I::FreeAssAlgIdeal{T}) where T
  return I.gens.Ox::parent_type(T)
end

function ngens(a::FreeAssAlgIdeal)
  return length(a.gens)
end

function gen(a::FreeAssAlgIdeal{T}, i::Int) where T
  return a.gens[Val(:O), i]
end

function gens(a::FreeAssAlgIdeal{T}) where T
  return T[gen(a,i) for i in 1:ngens(a)]
end


function Base.:+(a::FreeAssAlgIdeal{T}, b::FreeAssAlgIdeal{T}) where T
  R = base_ring(a)
  @assert R == base_ring(b) "parent mismatch"
  return ideal(R, vcat(gens(a), gens(b)))
end

function Base.:*(a::FreeAssAlgIdeal{T}, b::FreeAssAlgIdeal{T}) where T
  R = base_ring(a)
  @assert R == base_ring(b) "parent mismatch"
  return ideal(R, [i*j for i in gens(a) for j in gens(b)])
end

@doc raw"""
    ideal_membership(a::FreeAssAlgElem, I::FreeAssAlgIdeal, deg_bound::Int)

Return `true` if calculations with intermediate degrees bounded by `deg_bound`
prove that $a$ is in $I$. Otherwise, a return of `false` indicates an
inconclusive answer, but larger `deg_bound`s give more confidence in a negative
answer.
"""
function ideal_membership(a::FreeAssAlgElem, I::FreeAssAlgIdeal, deg_bound::Int)
  R = parent(a)
  @assert R == base_ring(I)
  groebner_assure(I, deg_bound)
  singular_assure(I.gb, deg_bound)
  Sx = base_ring(I.gb.S)
  return Singular.iszero(reduce(Sx(a), I.gb.S))
end

function Base.in(a::FreeAssAlgElem, I::FreeAssAlgIdeal, deg_bound::Int)
  return ideal_membership(a, I, deg_bound)
end

function (R::Singular.LPRing)(a::FreeAssAlgElem)
  B = MPolyBuildCtx(R)
  for (c, e) in zip(coefficients(a), exponent_words(a))
    push_term!(B, base_ring(R)(c), e)
  end
  return finish(B)
end

# ensure we have singular data with degree bound at least deg_bound
function singular_assure(I::IdealGens{T}, deg_bound::Int) where T <: FreeAssAlgElem
  deg_bound = max(deg_bound, 1)
  if !isdefined(I.gens, :S) || (isdefined(I.gens, :S) && I.Sx.deg_bound < deg_bound)
    oscar_assure(I)
    if !isdefined(I.gens, :S)
      # first time around: make sure the polys fit in the singular ring
      deg_bound = max(deg_bound, mapreduce(total_degree, max, I.O, init = deg_bound))
    end
    I.Sx = Singular.FreeAlgebra(singular_coeff_ring(coefficient_ring(I.Ox)),
                                symbols(I.Ox),
                                deg_bound)[1]
    I.S = Singular.Ideal(I.Sx, elem_type(I.Sx)[I.Sx(x) for x in I.O])
  end
  if I.isGB
    I.S.isGB = true
  end
end

# ensure we have singular groebner data with degree bound at least deg_bound
function groebner_assure(I::FreeAssAlgIdeal, deg_bound::Int)
  if !isdefined(I, :gb) || (isdefined(I.gb.gens, :Sx) && I.gb.Sx.deg_bound < deg_bound)
    singular_assure(I.gens, deg_bound)
    I.gb = IdealGens(I.gens.Ox, Singular.std(I.gens.S))
    I.gb.isGB  = true
  end
  return I.gb
end


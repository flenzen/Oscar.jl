# Detinko, Flannery, O'Brien "Recognizing finite matrix groups over infinite
# fields", Section 4.2
function _isomorphic_group_over_finite_field(matrices::Vector{T}) where T <: MatrixElem{<: Union{fmpq, nf_elem}}
   @assert !isempty(matrices)

   # One should probably check whether all matrices are n by n (and invertible
   # and such ...)

   K = base_ring(matrices[1])
   n = nrows(matrices[1])

   Fq, matrices_Fq, OtoFq = good_reduction(matrices, 2)

   G = MatrixGroup(n, Fq, matrices_Fq)
   N = order(G)
   if !isdivisible_by(Hecke._minkowski_multiple(K, n), N)
      error("Group is not finite")
   end

   G_to_fin_pres = GAP.Globals.IsomorphismFpGroupByGenerators(G.X, GapObj([ g.X for g in gens(G) ]))
   F = GAP.Globals.Range(G_to_fin_pres)
   rels = GAP.Globals.RelatorsOfFpGroup(F)

   gens_and_invsF = [ g for g in GAP.Globals.FreeGeneratorsOfFpGroup(F) ]
   append!(gens_and_invsF, [ inv(g) for g in GAP.Globals.FreeGeneratorsOfFpGroup(F) ])
   matrices_and_invs = copy(matrices)
   append!(matrices_and_invs, [ inv(M) for M in matrices ])
   for i = 1:length(rels)
      M = GAP.Globals.MappedWord(rels[i], GapObj(gens_and_invsF), GapObj(matrices_and_invs))
      if !isone(M)
         error("Group is not finite")
      end
   end
   return G, G_to_fin_pres, F, OtoFq
end

function good_reduction(matrices::Vector{T}, p::Int = 2) where T <: MatrixElem{<: Union{fmpq, nf_elem}}
   while true
      p = next_prime(p)
      b, Fq, matrices_Fq, OtoFq = test_modulus(matrices, p)
      b && return Fq, matrices_Fq, OtoFq
   end
end

# Small helper function to make the reduction call uniform
function _reduce(M::MatrixElem{nf_elem}, OtoFq::Hecke.NfOrdToFqMor)
  O = domain(OtoFq)
  Fq = codomain(OtoFq)
  return matrix(Fq, [ OtoFq(O(numerator(a)))//OtoFq(O(denominator(a))) for a in M])
end

function _reduce(M::MatrixElem{fmpq}, Fp)
  return map_entries(Fp, M)
end

function isomorphic_group_over_finite_field(G::MatrixGroup{T}) where T <: Union{fmpq, nf_elem}
   matrices = map(x -> x.elm, gens(G))

   Gp, GptoF, F, OtoFq = _isomorphic_group_over_finite_field(matrices)

   img = function(x)
     return Gp(_reduce(x.elm, OtoFq))
   end

   gen = gens(G)

   preimg = function(y)
     return GAP.Globals.MappedWord(GAP.Globals.UnderlyingElement(GAP.Globals.Image(GptoF, map_entries(Gp.ring_iso, y.elm))),
                                   GAP.Globals.FreeGeneratorsOfFpGroup(F),
                                   GAP.GapObj(gen))
   end

   return Gp, MapFromFunc(img, preimg, G, Gp)
end

# Detinko, Flannery, O'Brien "Recognizing finite matrix  groups over infinite
# fields", Section 3.1 claims that any prime != 2 not dividing any denominator
# of the matrices and their inverses (!) works, i.e. the projection is either
# an isomorphism or, if it is not injective, then the group generated by
# matrices cannot be finite.
function test_modulus(matrices::Vector{T}, p::Int) where T <: MatrixElem{fmpq}
   Fp = GF(p)
   matrices_Fp = Vector{AbstractAlgebra.MatElem{elem_type(Fp)}}(undef, length(matrices))
   if p == 2
      return false, Fp, matrices_Fp, Fp
   end

   for M in matrices
      for i = 1:nrows(M)
         for j = 1:ncols(M)
            if iszero(M[i, j])
               continue
            end

            if mod(denominator(M[i, j]), p) == 0
               return false, Fp, matrices_Fp, Fp
            end
         end
      end
   end
   # I don't want to invert everything in char 0, so I just check whether the
   # matrices are still invertible mod p.
   for i = 1:length(matrices)
      matrices_Fp[i] = map_entries(Fp, matrices[i])
      if rank(matrices_Fp[i]) != nrows(matrices_Fp[i])
         return false, Fp, matrices_Fp, Fp
      end
   end

   return true, Fp, matrices_Fp, Fp
end

# Detinko, Flannery, O'Brien "Recognizing finite matrix  groups over infinite
# fields", Section 3.2 claims that any prime != 2 not dividing the discriminant
# of the defining polynomial and not dividing any denominator of the matrices
# and their inverses (!) works, i.e. the projection is either
# an isomorphism or, if it is not injective, then the group generated by
# matrices cannot be finite.
function test_modulus(matrices::Vector{T}, p::Int) where T <: MatrixElem{nf_elem}
   @assert length(matrices) != 0
   K = base_ring(matrices[1])
   matrices_Fq = Vector{fq_mat}(undef, length(matrices))
   if p == 2
      return false, FiniteField(fmpz(p), 1, "a")[1], matrices_Fq, Hecke.NfOrdToFqMor()
   end
   O = EquationOrder(K)
   if mod(discriminant(O), p) == 0
      return false, FiniteField(fmpz(p), 1, "a")[1], matrices_Fq, Hecke.NfOrdToFqMor()
   end
   for M in matrices
      for i = 1:nrows(M)
         for j = 1:ncols(M)
            if iszero(M[i, j])
               continue
            end

            if mod(denominator(M[i, j]), p) == 0
               return false, FiniteField(fmpz(p), 1, "a")[1], matrices_Fq, Hecke.NfOrdToFqMor()
            end
         end
      end
   end

   # p is does not divide disc(O), so it's not an index divisor, so we don't
   # have to work in the maximal order here.
   P = prime_ideals_over(O, p)
   Fq, OtoFq = ResidueField(O, P[1])
   # I don't want to invert everything in char 0, so I just check whether the
   # matrices are still invertible mod p.
   for i = 1:length(matrices)
      matrices_Fq[i] = matrix(Fq, [ OtoFq(O(numerator(a)))//OtoFq(O(denominator(a))) for a in matrices[i] ])
      if rank(matrices_Fq[i]) != nrows(matrices_Fq[i])
         return false, Fq, matrices_Fq, Hecke.NfOrdToFqMor()
      end
   end

   return true, Fq, matrices_Fq, OtoFq
end

# Returns the largest possible order of a finite subgroup of GL(n, QQ) (equivalently:
# of GL(n, ZZ), as any finite subgroup of GL(n, QQ) is conjugate to a subgroup of GL(n, ZZ)).
# Always return a fmpz, only the orders for n <= 16 would fit into an Int64.
#
# This relies on results in a preprint "The orders of finite linear groups" by
# W. Feit (1995), possibly published as mathscinet.ams.org/mathscinet-getitem?mr=1484185
# in the "Proceedings of the First Jamaican Conference on Group Theory and its
# Applications", 1996. However, it seems basically impossible to that paper.
# Geoff Robinson claims to have the preprint and posted the relevant information
# at mathoverflow.net/questions/168292/maximal-order-of-finite-subgroups-of-gln-z .
# The table is also repeated in [BDEPS04] where the authors however state that
# Feit does not actually provide a proof, and in any case relies heavily on unpublished
# work by Weisfeiler.
#
# Since this leaves the result in doubt, we do not currently actually use this code,
# and instead resort to the Minkowski bound and its generalization. Luckily for us,
# Hecke already implements these.
const max_ords = [ 2, 12, 48, 1152, 3840, 103680, 2903040, 696729600, 1393459200, 8360755200 ]
function largest_order_of_finite_linear_group_qq(n::Int)
   @assert n >= 0
   n <= 10 && return fmpz(max_ords[n])
   # For n > 10, we can use 2^n*n!
   return factorial(fmpz(n)) << n
end

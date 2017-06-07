include("matrixops.jl")
using MatrixOps

module Motif

"""
`symmetric_stochastic_block_model(m,k,p,q) -> SparseMatrixCSC, sets`
"""
function symmetric_stochastic_block_model(m,k,p,q)
  n = m*k
  A = sprand(n,n,q)
  offset = 1
  sets = Vector{Vector{Int}}()
  for i=1:k
    A[offset:(offset+m-1), offset:(offset+m-1)] = sprand(m,m,p)
    push!(sets, offset:(offset+m-1))
    offset += m
  end
  A = spones(triu(A,1))
  A = A + A';
  return A, sets
end

"""
ssbm_round(M, k, p, q) -> SparseMatrixCSC
M is a dense matrix of probabilities
"""
function ssbm_round(M,m,k,p,q)
  bits = (M .<= q) # get the raw one
  offset = 1
  for i=1:k
    bits[offset:(offset+m-1), offset:(offset+m-1)] = M[offset:(offset+m-1), offset:(offset+m-1)] .<= p
    offset += m
  end
  A = spones(triu(sparse(bits),1))
  return A + A'
end

"""
  `score_set(iterable,m,k) -> Float`

return the best accuracy for any of the k sets of size m given by the standard indices.
"""
function score_set(iter,m,k)
  best = 0.0
  offset = 1
  for i=1:k
    acc = length(intersect(iter, offset:(offset+m-1)))/m
    if acc >= best
      best = acc
    end
    offset += m
  end
  return best
end

using MatrixOps
using MatrixNetworks

function network_matrices(A::SparseMatrixCSC)
  n = size(A,1)
  mats = Dict{Symbol,SparseMatrixCSC}()
  ops = Dict{Symbol,MatrixCFcn}()
  mats[:Adjacency] = A

  d = sum(A,2)
  dhalf = sqrt(d)
  dhalfnorm = dhalf/norm(dhalf)
  K = diagm(sparse(d))-A
  mats[:Kirchoff] = K
  mats[:CombinatorialLaplacian] = K

  ei,ej,ev = findnz(A)
  W = sparse(ei,ej,ev./(dhalf[ei].*dhalf[ej]),n,n);
  mats[:WeightedAdjacency] = W

  mats[:NormalizedLaplacian] = speye(n) - W

  #mats[:ShiftedLaplacian] = speye(n) + W
  ops[:ShiftedLaplacian] = MatrixCFcn{eltype(W)}(n, n,
    (output, x) -> begin
      A_mul_B!(output, W, x)
      output[:] += -2*(dhalfnorm'*x)[1]*dhalfnorm
      output[:] += x
    end,
    (output, x) -> begin
      A_mul_B!(output, W, x)
      output[:] += -2*(dhalfnorm'*x)[1]*dhalfnorm
      output[:] += x
    end
  )

  return mats, ops
end

end # end module

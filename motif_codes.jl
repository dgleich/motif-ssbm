include("matrixops.jl")
using MatrixOps

module Motif

"""
`symmetric_stochastic_block_model(m,k,p,q) -> SparseMatrixCSC, sets`
"""
function symmetric_stochastic_block_model(m,k,p,q)
  n = m*k
  A = sprand(Bool,n,n,q)
  offset = 1
  sets = Vector{Vector{Int}}()
  for i=1:k
    A[offset:(offset+m-1), offset:(offset+m-1)] = sprand(Bool,m,m,p)
    push!(sets, offset:(offset+m-1))
    offset += m
  end
  A = spones(dropzeros!(triu(A,1)))
  A = float(A + A')
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
  best_set_and_score(iter,m,k)[2]
end

function best_set_and_score(iter,m,k)
  best = 0.0
  bestk = 0
  offset = 1
  for i=1:k
    acc = length(intersect(iter, offset:(offset+m-1)))/m
    if acc >= best
      best = acc
      bestk = i
    end
    offset += m
  end
  return bestk,best
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

function clique4_weighted{T}(A::SparseMatrixCSC{T,Int64})
    C = max(A, A')
    order = sortperm(vec(sum(C, 1)))
    B = tril(C[order, order])
    C = deepcopy(A)[order, order]
    n = size(B, 1)

    x = Int64[];
    y = Int64[];
    z = Int64[];
    for i = 1:n
        nbrs = find(B[:,i])
        nnbr = length(nbrs)
        for jj = 1:nnbr
            for kk = (jj+1):nnbr
                for ll = (kk+1):nnbr
                    j = nbrs[jj]
                    k = nbrs[kk]
                    l = nbrs[ll]
                    if i == j || i == k || i == l || j == k  || j == l || k == l
                        continue
                    end
                    if C[j, k] != 0 && C[j, l] != 0 && C[k, l] != 0
                        push!(x, i, i, i, j, j, k)
                        push!(y, j, k, l, k, l, l)
                        push!(z, 1, 1, 1, 1, 1, 1)
                    end
                end
            end
        end
        if (i % 100000 == 0)
            @show i
        end
    end

    # Make sure it is square
    maxdim = maximum(size(A))
    push!(x, maxdim)
    push!(y, maxdim)
    push!(z, 0)
    W = sparse(x, y, z)

    # Reorder
    rev_order = zeros(Int64, n)
    rev_order[order] = 1:n
    W = W[rev_order, rev_order]
    return W + W'
end


end # end module

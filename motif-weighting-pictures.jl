## Setup
include("motif_codes.jl")
using MatrixNetworks
using Plots


## Animate the ssbm

function my_symmetric_stochastic_block_model(m,k,p,q)
  n = m*k
  A = sprand(Bool,n,n,q)
  offset = 1
  sets = Vector{Vector{Int}}()
  @show p, q
  for i=1:k
    A[offset:(offset+m-1), offset:(offset+m-1)] = sprand(Bool,m,m,p)
    push!(sets, offset:(offset+m-1))
    offset += m
  end
  A = spones(dropzeros!(triu(A,1)))
  A = float(A + A')
  return A, sets
end

m = 200
k = 5
p = 0.3

srand(1)
M = rand(m*k,m*k)
npts = 50

pyplot(size=(1600,800))
μs = [linspace(0,0.6,npts);0.6*ones(15);linspace(0.6,0.9,npts);0.9*ones(20)]

q = 0.13

srand(1)
A,sets = my_symmetric_stochastic_block_model(m,k,p,q)
M = (A*A).*A

function expand_nonzeros(M)
  ei,ej,ev = findnz(M)
  is = zeros(Int,0)
  js = zeros(Int,0)
  for nzi in eachindex(ei)
      for r = 1:Int(ev[nzi])
        push!(is, ei[nzi]-1)
        push!(js, ej[nzi]-1)
      end
  end
  return is,js
end

plot(
  #heatmap(full(A),yflip=true,leg=false,border=false,xticks=[],yticks=[]),
  histogram2d(map( (i,j) -> (i-1,j-1), findnz(A)[1:2]...), nbins=50, yflip=true,leg=false,border=false,xticks=[],yticks=[]),
#heatmap(full(((A*A).*A)).^2,yflip=true,leg=false,border=false,xticks=[],yticks=[])
  #heatmap(full(M.*(M .>= 1)),yflip=true,leg=false,border=false,xticks=[],yticks=[])
  histogram2d(expand_nonzeros(M)..., nbins=50, yflip=true,leg=false,border=false,xticks=[],yticks=[]),
)
gui()

##
pyplot(size=(800,800))
histogram2d(map( (i,j) -> (i-1,j-1), findnz(A)[1:2]...), nbins=50, yflip=true,leg=false,border=false,xticks=[],yticks=[])
savefig("ssbm-$m-$k-$p-$q-edge.png")

histogram2d(expand_nonzeros(M)..., nbins=50, yflip=true,leg=false,border=false,xticks=[],yticks=[])
savefig("ssbm-$m-$k-$p-$q-motif.png")

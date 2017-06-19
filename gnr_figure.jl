##
using Plots
using PlotRecipes
using NearestNeighbors
using Roots

##
function gnr(n,r)
  xy = rand(2,n)
  T = BallTree(xy)
  idxs = inrange(T, xy, r)
  # form the edges for sparse
  ei = Int[]
  ej = Int[]
  for i=1:n
    for j=idxs[i]
      if i > j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  return xy, ei, ej
end

srand(2)
xy, ei, ej = gnr(50,0.23)

pyplot(size=(300,300))
graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
  markercolor=:black, markerstrokecolor=:white,
  markersize=4, linecolor=1, linealpha=0.8, linewidth=0.7,
  axis_buffer=0.02, background=nothing)
gui()
savefig("gnr-figure.pdf")

##
A = sparse(ei,ej,1.0,50,50)
A = A+A'
spy(A,markersize=5,background=nothing,border=:black,xticks=[],yticks=[])
xlims!(0,51)
ylims!(0,51)
gui()
#heatmap(full(A),markersize=5,background=nothing,border=:black,xticks=[],yticks=[])
savefig("gnr-figure-adj.pdf")

##
function gnk(n,k)
  xy = rand(2,n)
  T = BallTree(xy)
  idxs = knn(T, xy, k)[1]
  # form the edges for sparse
  ei = Int[]
  ej = Int[]
  for i=1:n
    for j=idxs[i]
      if i != j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  return xy, ei, ej
end

srand(1)
pyplot(size=(300,300))
xy, ei, ej = gnk(50,4)
A = sparse(ei,ej,1.0,50,50)
B = A.*A'
U = A - B

graphplot(findnz(triu(B,1))[1:2]..., x = xy[1,:], y = xy[2,:],
  markercolor=:black, markerstrokecolor=:white, curves=false,
  markersize=0, linecolor=1, linealpha=0.8, linewidth=0.7,
  axis_buffer=0.02, background=nothing)

graphplot!(findnz(U)[1:2]..., x = xy[1,:], y = xy[2,:],
  markercolor=:black, markerstrokecolor=:white, curves=true, curvature_scalar = -0.4,
  markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7, arrow = 0.3, shorten=0.85,
  axis_buffer=0.02, background=nothing)

gui()
savefig("gnk-figure-adj.pdf")

##
A = sparse(ei,ej,1.0,50,50)
spy(A,markersize=5,background=nothing,border=:black,xticks=[],yticks=[])
xlims!(0,51)
ylims!(0,51)
gui()
#heatmap(full(A),markersize=5,background=nothing,border=:black,xticks=[],yticks=[])
savefig("gnk-figure-adj.pdf")

##
A = sparse(ei,ej,1.0,50,50)
B = A.*A'
U = dropzeros(A - B)
scatter(findnz(B)[1:2]...,markerstrokewidth=0, color=:black,
  markersize=5,background=nothing,border=:black,xticks=[],yticks=[],leg=false)
scatter!(findnz(U)[1:2]...,markerstrokewidth=0,color=:firebrick,
    markersize=5,background=nothing,border=:black,xticks=[],yticks=[],leg=false)
#spy(B,markersize=5,background=nothing,border=:black,xticks=[],yticks=[],marker_z=:red)
#spy!(U,markersize=5,color=2)
xlims!(0,51)
ylims!(0,51)
gui()
savefig("gnk-figure-adj-colors.pdf")

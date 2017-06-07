## Setup
include("motif_codes.jl")
using MatrixNetworks
using Plots


## Animate the ssbm

# Version 1
m = 50
k = 10
p = 0.20

srand(1)
M = rand(m*k,m*k)
npts = 50

pyplot(size=(400,400))
μs = [linspace(0,0.6,npts);0.6*ones(15);linspace(0.6,0.9,npts);0.9*ones(20)]
anim = @animate for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))

  A = Motif.ssbm_round(M,m,k,p,q)
  MA = (A*A).*A
  #plt = spy(MA, layout=l, title="μ = $(μ)", legend=false)
  ei,ej,vals = findnz(A)
  histogram2d( ei-1,ej-1,vals, yflip=true,nbins=100, legend=false, title=@sprintf("μ = %.3f", μ))
  xticks!(0:m:m*k)
  yticks!(0:m:m*k)


  #spy!(plt[3], A)
end
gif(anim, "ssbm-adjacency.gif", fps=15)


##

srand(1)
M = rand(m*k,m*k)
npts = 50

pyplot(size=(400,400))
μs = [linspace(0,0.6,npts);0.6*ones(15);linspace(0.6,0.9,npts);0.9*ones(20)]
anim = @animate for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))

  A = Motif.ssbm_round(M,m,k,p,q)
  MA = (A*A).*A
  #plt = spy(MA, layout=l, title="μ = $(μ)", legend=false)
  ei,ej,vals = findnz(MA)
  histogram2d( ei-1,ej-1,vals, yflip=true,nbins=100, legend=false, title=@sprintf("μ = %.3f", μ))
  xticks!(0:m:m*k)
  yticks!(0:m:m*k)


  #spy!(plt[3], A)
end
gif(anim, "ssbm-triangle-motif-adjacency.gif", fps=15)

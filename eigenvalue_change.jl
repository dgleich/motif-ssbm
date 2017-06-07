## This was just some initial work as I was trying to figure out the best presentation
# of this material. Follow up codes with more insightful things.
# * powermethod_steps.jl
# * powermethod_animation.jl


## Setup
include("motif_codes.jl")
using MatrixNetworks
using Plots

## How do the eigenvalues change
m = 50
k = 10
p = 0.40
μ = 0.6
q = μ*p/((1-μ)*(k-1))

srand(1)
A,grps = Motif.symmetric_stochastic_block_model(m,k,p,q)
M = (A*A).*A
@show is_connected(M)
@assert (nnz(M-M')==0)
n = size(A,1)
##
Amats = Motif.network_matrices((A))[1]
Mmats = Motif.network_matrices((M))[1]
AAlams,AAvecs = eig(full(Amats[:Adjacency]))
AMlams,AMvecs = eig(full(Mmats[:Adjacency]))

scatter(1:n, AAlams, label="Edge")
scatter!(1:n, AMlams, label="Triangle")
gui()

##
#scatter(1:n, [ AAvecs[:,end-k+1:end] AMvecs[:,end-k+1:end]],layout=2)
Plots.PyPlot.plt[:close]("all")
scatter(1:n, [ AAvecs[:,end-k+1:end] ])
gui()
Plots.PyPlot.plt[:figure]()
scatter(1:n, [ AMvecs[:,end-k+1:end] ])
gui()

##
LAlams,LAvecs = eig(full(Amats[:NormalizedLaplacian]))
LMlams,LMvecs = eig(full(Mmats[:NormalizedLaplacian]))

scatter(1:n, LAlams, label="Edge")
scatter!(1:n, LMlams, label="Triangle")
gui()

## Show the Fiedler vectors
scatter(1:n, [fiedler_vector(Amats[:Adjacency])[1] fiedler_vector(Mmats[:Adjacency])[1]])
gui()
## Show the smallest few eigenvalues
#scatter(1:n, [ AAvecs[:,end-k+1:end] AMvecs[:,end-k+1:end]],layout=2)
Plots.PyPlot.plt[:close]("all")
scatter(1:n, [ diagm(1./vec(sqrt(sum(Amats[:Adjacency],2))))*LAvecs[:,1:k] ])
gui()
Plots.PyPlot.plt[:figure]()
scatter(1:n, [ diagm(1./vec(sqrt(sum(Mmats[:Adjacency],2))))*LMvecs[:,1:k] ])
gui()

## This experiment suggests it'd be cool to animate SSBM as we vary mu, let's try that!
srand(1)
m = 50
k = 10
p = 0.40

M = rand(m*k,m*k)
μ = 0.6
A = Motif.ssbm_round(M,m,k,p,q)
npts = 10

μs = linspace(0.1,0.8,npts)
for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))
  A = Motif.ssbm_round(M,m,k,p,q)
  spy(A; title="μ = $(μ)", legend=false)
  gui()
  sleep(0.1)
end

##
srand(1)
m = 50
k = 10
p = 0.40

M = rand(m*k,m*k)
μ = 0.6
A = Motif.ssbm_round(M,m,k,p,q)
npts = 10

μs = linspace(0.1,0.8,npts)
for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))
  A = Motif.ssbm_round(M,m,k,p,q)
  MA = (A*A).*A
  spy(MA; title="μ = $(μ)", legend=false)
  gui()
  sleep(0.1)
end

## Let's animate the Fiedler vector
srand(1)
m = 50
k = 10
p = 0.40

M = rand(m*k,m*k)
μ = 0.6
A = Motif.ssbm_round(M,m,k,p,q)
npts = 10

μs = linspace(0.1,0.8,npts)
for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))
  A = Motif.ssbm_round(M,m,k,p,q)
  l = @layout([a{0.9h}; b])
  plt = spy(A, layout=l, title="μ = $(μ)", legend=false)
  scatter!(plt[2],1:n,fiedler_vector(A)[1])
  gui()
  sleep(0.1)
end


## Let's animate the Fiedler vector
srand(2)
m = 50
k = 10
p = 0.40

M = rand(m*k,m*k)
μ = 0.6
A = Motif.ssbm_round(M,m,k,p,q)
npts = 25
accs = zeros(npts,2)

μs = linspace(0.1,0.8,npts)
for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))
  A = Motif.ssbm_round(M,m,k,p,q)
  MA = (A*A).*A
  l = @layout([a{0.9h} c; b d; e{0.1h}])
  #plt = spy(MA, layout=l, title="μ = $(μ)", legend=false)
  plt = plot(layout=l, legend=false)
  ei,ej,vals = findnz(MA)
  histogram2d!(plt[1], ei-1,ej-1,vals, yflip=true,nbins=100)
  z = fiedler_vector(MA)[1]
  z *= sign(z[indmax(abs(z))])
  perm = sortperm(z)
  acc_motif = max(Motif.score_set(perm[1:m],m,k), Motif.score_set(perm[end-m+1:end],m,k))
  scatter!(plt[3],1:n,z, legend=false)

  ei,ej,vals = findnz(A)
  histogram2d!(plt[2], ei-1,ej-1,vals, yflip=true,nbins=100)
  z = fiedler_vector(A)[1]
  z *= sign(z[indmax(abs(z))])
  perm = sortperm(z)
  acc_edge = max(Motif.score_set(perm[1:m],m,k), Motif.score_set(perm[end-m+1:end],m,k))
  scatter!(plt[4],1:n,z, legend=false)

  accs[i,1] = acc_motif
  accs[i,2] = acc_edge

  plot!(plt[5],μs, accs[:,1], legend=true,label="motif")
  plot!(plt[5],μs, accs[:,2], legend=true,label="edge")

  #spy!(plt[3], A)
  gui()
  sleep(0.1)
end

## Let's

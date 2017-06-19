## Setup
include("motif_codes.jl")
using MatrixNetworks
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false" # turn off Atom plotpane
using Plots
pyplot(size=(800,600))
using StatsBase

## We want to look at
# new types of distrubtions that fit better.
# Let's start with just ER and look at the eigenvalues of the motifs
# Okay, this has awful scaling properties, sigh.
# Back to Laplacians
begin
  p = 0.15
  q = 0.1
  m = 150
  k = 1
  ntris = zeros(0)
  lams = zeros(0)
  lams_A = zeros(0)
  for t=1:100
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A
    push!(ntris, sum(M))
    push!(lams, eigvals!(full(M))...)
    push!(lams_A, eigvals!(full(A))...)
  end

  h = fit(Histogram, lams; nbins = 100)
  h = normalize(h)
  plot(h.edges, h.weights)

  h = fit(Histogram, lams_A; nbins = 100)
  h = normalize(h)
  plot!(h.edges, h.weights)

  h = fit(Histogram, lams_A.^2; nbins = 100)
  h = normalize(h)
  plot!(h.edges, h.weights)

  h = fit(Histogram, rand(lams_A.^3,10000); nbins = 100)
  h = normalize(h)
  plot!(h.edges, h.weights)

  etris = (mean(ntris)./(m*k))
  gui()
end


##
begin
  p = 0.2
  q = 0.0
  m = 50
  k = 12
  ntris = zeros(0)
  lams = zeros(0)
  lams_A = zeros(0)
  lams_A2 = zeros(0)
  for t=1:50
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A
    push!(ntris, sum(M))
    push!(lams, eigvals!(full(Motif.network_matrices(M)[1][:NormalizedLaplacian]))...)
    push!(lams_A, eigvals!(full(Motif.network_matrices(A)[1][:NormalizedLaplacian]))...)
    push!(lams_A2, eigvals!(full(Motif.network_matrices(A*A)[1][:NormalizedLaplacian]))...)
  end
  mprho = 3.0/(mean(ntris)./(m*k))
end
##

begin
  h = fit(Histogram, (abs(lams)); nbins = 100)
  h = normalize(h)
  plot(h.edges, h.weights)

  h = fit(Histogram, (abs(lams)).^2; nbins = 100)
  h = normalize(h)
  plot!(h.edges, h.weights)

  #h = fit(Histogram, lams_A; nbins = 100)
  #h = normalize(h)
  #plot!(h.edges, h.weights)


  # Insight, this seems to be a good example, in particular, the
  # lower bound seems to be fairly good.
  pbar = (q*(k)*(k-1)+p*k)/(k^2)
  cla = (1-sqrt(3/(m*k*pbar)))*(1-sqrt(mprho))^2
  clb = (1+sqrt(3/(m*k*pbar)))*(1+sqrt(mprho))^2
  vline!([cla,clb])

  h = fit(Histogram, 2.0-lams_A.^2+sqrt(mprho)/2; nbins = 100)
  h = normalize(h)
  plot!(h.edges, h.weights)

  #h = fit(Histogram, lams_A2; nbins = 100)
  #h = normalize(h)
  #plot!(h.edges, h.weights)


  #@show mprho = 3.0/(mean(ntris)./(m*k))

  gui()
end

##
# We have the region 0.641 to 1.25
# for the eigenvalues of the motif matrix
# a good scaling measure seems to be 0.0137
mprho = 0.0137
@show 2-(1-sqrt(mprho))^2.5
@show 2-(1+sqrt(mprho))



##
#=
mprho = 5.0/(mean(ntris)./(m*k))
mpa = (1-sqrt(mprho))^2
mpb = (1+sqrt(mprho))^2
plot!(h.edges[1], map(lam -> begin
  lam = 2.0-lam
  if mpa <= lam <= mpb
    return sqrt((mpb - lam)*(lam-mpa)) /(2π*lam*mprho)
  else
    return 0.0
  end
end, h.edges[1]))
=#

#=
mprho = 0.01
mpa = (1-sqrt(mprho))^2
mpb = (1+sqrt(mprho))^2
@show 2.0-mpa, 2.0-mpb
plot!(collect(linspace(0,2,100)), map(lam -> begin
  lam = 2.0-lam
  if mpa <= lam <= mpb
    return sqrt((mpb - lam)*(lam-mpa)) /(2π*lam*mprho)
  else
    return 0.0
  end
end, collect(linspace(0,2,100))))
=#

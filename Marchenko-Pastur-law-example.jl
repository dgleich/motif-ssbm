## Setup
include("motif_codes.jl")
using MatrixNetworks
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false" # turn off Atom plotpane
using Plots
pyplot(size=(800,600))
using StatsBase

## We want to look at
begin
  p = 0.2
  q = 0.0
  m = 100
  k = 6
  ntris = zeros(0)
  lams = zeros(0)
  for t=1:200
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A
    push!(ntris, sum(M))
    push!(lams, eigvals!(full(Motif.network_matrices(M)[1][:NormalizedLaplacian]))...)
  end

  h = fit(Histogram, lams; nbins = 100)
  h = normalize(h)
  plot(h.edges, h.weights)
  mprho = 3.0/(mean(ntris)./(m*k))
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
  gui()
end

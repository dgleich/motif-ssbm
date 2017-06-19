# In this one, we are looking at various conditioning and sensitivity measures
# of the problem.

## Setup
include("motif_codes.jl")
using MatrixNetworks
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false" # turn off Atom plotpane
using Plots
using StatsBase
## Based on eigenvalue_sensitivity, we found that the Motif laplacian had a bigger
# gap deeper in the spectrum, and this seemed to explain the convergence of the
# power method. So we'd like to see this as a 2d histogram.

m = 50
k = 10
p = 0.2

srand(1)
npts = 100
μs = linspace(0.1,0.8,npts)

nhistbins = 50
histbins = linspace(0,2,nhistbins+1)
histpts = diff(histbins)+histbins[1:end-1]

lams_A = zeros(m*k, npts)
lams_M = zeros(m*k, npts)

hist_A = zeros(nhistbins, npts)
hist_M = zeros(nhistbins, npts)

for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))

  A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
  M = (A*A).*A

  lams_A[:,i] = eigvals!(full(Motif.network_matrices(A)[1][:NormalizedLaplacian]))
  lams_M[:,i] = eigvals!(full(Motif.network_matrices(M)[1][:NormalizedLaplacian]))

  hist_A[:,i] = fit(Histogram, lams_A[:,i, ], histbins).weights
  hist_M[:,i] = fit(Histogram, lams_M[:,i, ], histbins).weights

  # generate additional trials to smooth out the histograms
  for t=1:50
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A

    hist_A[:,i] += fit(Histogram, eigvals!(full(Motif.network_matrices(A)[1][:NormalizedLaplacian])), histbins).weights
    hist_M[:,i] += fit(Histogram, eigvals!(full(Motif.network_matrices(M)[1][:NormalizedLaplacian])), histbins).weights
  end
end

##
nevals = 10
pyplot(size=(300,300))
Mus = ones(nevals,1)*μs'
#plt = plot(layout=2, leg=false, ylim=(0.1,0.8), xlim=(0,2))
plt = plot(layout=1, leg=false, ylim=(0.1,0.8), xlim=(0,2))
contour!(plt[1],  histpts, μs, hist_A', leg=false, fill=true, levels=8),
scatter!(plt[1], vec(lams_A[2:2+nevals-1,:]), vec(Mus),
  markerstrokewidth=0, markeralpha=0.5, label="", leg=false, border=false),
xlabel!("Eigenvalue")
ylabel!("μ")
#contour!(plt[2], μs, histpts, hist_M, fill=true, levels=8),
#scatter!(plt[2], vec(Mus), vec(lams_M[2:2+nevals-1,:]), markerstrokewidth=0, markeralpha=0.5, label=""),

gui()
savefig("edge-eigenvalues.pdf")

##
nevals = 10
pyplot(size=(300,300))
Mus = ones(nevals,1)*μs'
#plt = plot(layout=2, leg=false, ylim=(0.1,0.8), xlim=(0,2))
plt = plot(layout=1, leg=false, ylim=(0.1,0.8), xlim=(0,2))
contour!(plt[1],  histpts, μs, hist_M', leg=false, fill=true, levels=8),
scatter!(plt[1], vec(lams_M[2:2+nevals-1,:]), vec(Mus),
  markerstrokewidth=0, markeralpha=0.5, label="", leg=false, border=false),
xlabel!("Eigenvalue")
ylabel!("μ")
#contour!(plt[2], μs, histpts, hist_M, fill=true, levels=8),
#scatter!(plt[2], vec(Mus), vec(lams_M[2:2+nevals-1,:]), markerstrokewidth=0, markeralpha=0.5, label=""),

gui()
savefig("motif-eigenvalues.pdf")
##
nevals = 10
Mus = ones(nevals,1)*μs'
plt = plot(layout=2, leg=false, xlim=(0.1,0.8), ylim=(0,2))
contour!(plt[1], μs, histpts, hist_A, leg=false, fill=true, levels=8),
scatter!(plt[1], vec(Mus), vec(lams_A[2:2+nevals-1,:]),
  markerstrokewidth=0, markeralpha=0.5, label="", leg=false, border=false),

contour!(plt[2], μs, histpts, hist_M, fill=true, levels=8),
scatter!(plt[2], vec(Mus), vec(lams_M[2:2+nevals-1,:]), markerstrokewidth=0, markeralpha=0.5, label=""),

gui()

##
Mus = ones(m*k,1)*μs'
#histogram2d(vec(Mus), vec(lams_A),nxbins=npts)
plot(
  scatter(vec(Mus), vec(lams_A)),
  scatter(vec(Mus), vec(lams_M)),
  markerstrokewidth=0, markeralpha=0.2)
gui()

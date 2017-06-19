##
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
npts = 500
μs = linspace(0.0,0.8,npts)

nhistbins = 50
histbins = linspace(0,2,nhistbins+1)
histpts = diff(histbins)+histbins[1:end-1]

lams_A = zeros(m*k, npts)
lams_M = zeros(m*k, npts)

hist_A = zeros(nhistbins, npts)
hist_M = zeros(nhistbins, npts)

ntris = zeros(npts)
phi1 = zeros(npts)
phihalf = zeros(npts)
phirand = zeros(npts)

for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))

  A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
  M = (A*A).*A

  ntris[i] = sum(M) # this gives 6 times the number of triangles
  phi1[i] = sum(M[m+1:end,1:m])/sum(M[:,1:m])
  volhalf = sum(M[:,1:(m*floor(Int,k/2))])

  # compute conductance
  myphi = S -> begin
    vol = sum(M[:,S])
    cut = vol - sum(M[S,S])
    return cut/min(ntris[i]-vol,vol)
  end


  phihalf[i] = sum(
                M[(m*floor(Int,k/2)+1):end,1:(m*floor(Int,k/2))])/
               min(volhalf,ntris[i]-volhalf)

  phirand[i] = mean( map( x -> myphi(randperm(m*k)[1:(m*div(k,2))]), 1:100 ))

  lams_A[:,i] = eigvals!(full(Motif.network_matrices(A)[1][:NormalizedLaplacian]))
  lams_M[:,i] = eigvals!(full(Motif.network_matrices(M)[1][:NormalizedLaplacian]))

  hist_A[:,i] = fit(Histogram, lams_A[:,i, ], histbins).weights
  hist_M[:,i] = fit(Histogram, lams_M[:,i, ], histbins).weights

  # generate additional trials to smooth out the histograms
  for t=1:0
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A

    hist_A[:,i] += fit(Histogram, eigvals!(full(Motif.network_matrices(A)[1][:NormalizedLaplacian])), histbins).weights
    hist_M[:,i] += fit(Histogram, eigvals!(full(Motif.network_matrices(M)[1][:NormalizedLaplacian])), histbins).weights
  end
end

##
nevals = k
pyplot(size=(300,300))
Mus = ones(nevals,1)*μs'
#plt = plot(layout=2, leg=false, ylim=(0.1,0.8), xlim=(0,2))
plt = plot(layout=1, leg=false, ylim=(0.0,0.8), xlim=(0,2))
contour!(plt[1],  histpts, μs, hist_A', leg=false, fill=true, levels=8),
scatter!(plt[1], vec(lams_A[2:2+nevals-1,:]), vec(Mus),
  markerstrokewidth=0, markeralpha=0.5, label="", leg=false, border=false),
xlabel!("Eigenvalue")
ylabel!("μ")


Plots.abline!(1.0,0.0)
#vline!([1-sqrt(3./((m*k*p))),1+sqrt(3./((m*k*p)))])
#vline!([1-sqrt(3./((m*p))),1+sqrt(3./((m*p)))])

# This lower-bound comes from
# https://arxiv.org/pdf/1607.07069.pdf
# Theorem 23.3.7, using the average ER probability over the SBM.
lbs = map( μ -> begin
                  q = μ*p/((1-μ)*(k-1))
                  pbar = (q*(k)*(k-1)+p*k)/(k^2)
                  return 1-sqrt(3/(m*k*pbar))
                end,
           μs )
plot!(lbs, collect(μs))

#Plots.abline!(1.0,+0.05)
#contour!(plt[2], μs, histpts, hist_M, fill=true, levels=8),
#scatter!(plt[2], vec(Mus), vec(lams_M[2:2+nevals-1,:]), markerstrokewidth=0, markeralpha=0.5, label=""),

gui()

##
nevals = k
Mus = ones(nevals,1)*μs'
#plt = plot(layout=2, leg=false, ylim=(0.1,0.8), xlim=(0,2))
plt = plot(layout=1, leg=false, ylim=(0.0,0.8), xlim=(0,2))
contour!(plt[1],  histpts, μs, hist_M', leg=false, fill=true, levels=8),
scatter!(plt[1], vec(lams_M[2:2+nevals-1,:]), vec(Mus),
  markerstrokewidth=0, markeralpha=0.5, label="", leg=false, border=false),
xlabel!("Eigenvalue")
ylabel!("μ")
#contour!(plt[2], μs, histpts, hist_M, fill=true, levels=8),
#scatter!(plt[2], vec(Mus), vec(lams_M[2:2+nevals-1,:]), markerstrokewidth=0, markeralpha=0.5, label=""),

phis = map(
  μ -> begin
    q = μ*p/((1-μ)*(k-1))
    cut1 = m*(m*(k-1))*(m*(k-2))*q^3
    cut2 = m*(m*(k-1))*(m-1)*q^2*p
    return (cut1+cut2)/(m*(m-1)*(m-2)*p^3 + cut1+cut2)
  end, μs )
plot!(phis, collect(μs))

# Austin worked out this one.
phisa = map(
  μ -> begin
    q = μ*p/((1-μ)*(k-1))
    cut1 = p*q^2*binomial(m,2)*(k*m-1)
    cut2 = p*q^2*m*(k-1)*binomial(m,2)
    cut3 = q^3*m^3*binomial(k-1,2)
    vol1 = 3*p^3*binomial(m,3)
    vol2 = 2*cut1
    vol3 = cut2
    vol4 = cut3
    return (cut1+cut2+cut3)/(vol1+vol2+vol3+vol4)
  end, μs )
plot!(phisa, collect(μs))

lbs = map(
  μ -> begin
    q = μ*p/((1-μ)*(k-1))
    type1 = p^3*k*m*(m-1)*(m-2)
    type2 = 3*p*q^2*k*(m)*(m-1)*(k-1)*m
    type3 = 2*q^3*(k*m)*((k-1)*m)*((k-2)*m)
    expected_tris = ((type1+type2+type3))/(k*m)
    return 1-sqrt(3/(expected_tris))
  end, μs)
#plot!(lbs, collect(μs))
#plot!(2.0-lbs, collect(μs))


lbs3 = 1.0-sqrt(3./(ntris./(m*k)))

lbs_mp_cl = map( μ -> begin
    q = μ*p/((1-μ)*(k-1))
    type1 = p^3*k*m*(m-1)*(m-2)
    type2 = 3*p*q^2*k*(m)*(m-1)*(k-1)*m
    type3 = 2*q^3*(k*m)*((k-1)*m)*((k-2)*m)
    pbar = (q*(k)*(k-1)+p*k)/(k^2)
    expected_tris = ((type1+type2+type3))/(k*m)
    mprho = 3.0/(expected_tris)
    #return (1-sqrt(3.0/(m*k*pbar)))*(1-sqrt(mprho))^2
    #return (1-sqrt(mprho))^2
    return (1-sqrt(3.0/(m*k*pbar)))*(1-sqrt(mprho))^2
  end, μs)


lbs_mp = map( μ -> begin
    q = μ*p/((1-μ)*(k-1))
    type1 = p^3*k*m*(m-1)*(m-2)
    type2 = 3*p*q^2*k*(m)*(m-1)*(k-1)*m
    type3 = 2*q^3*(k*m)*((k-1)*m)*((k-2)*m)
    expected_tris = ((type1+type2+type3))/(k*m)
    mprho = 2.0/(expected_tris)
    mpa = 2-(1-sqrt(mprho))^2
    mpb = 2-(1+sqrt(mprho))^2
    return mpb
  end, μs)


ubs_mp = map( μ -> begin
      q = μ*p/((1-μ)*(k-1))
      type1 = p^3*k*m*(m-1)*(m-2)
      type2 = 3*p*q^2*k*(m)*(m-1)*(k-1)*m
      type3 = 2*q^3*(k*m)*((k-1)*m)*((k-2)*m)
      expected_tris = ((type1+type2+type3))/(k*m)
      mprho = 2.0/(expected_tris)
      mpa = 2-(1-sqrt(mprho))^2
      mpb = 2-(1+sqrt(mprho))^2
      return mpa
    end, μs)


#plot!((lbs+phisa)/2, collect(μs))
#plot!(lbs2, collect(μs))
#plot!(lbs3, collect(μs))
plot!(lbs_mp_cl, collect(μs))
plot!(ubs_mp, collect(μs))
#plot!(phi1, collect(μs))
#plot!(phihalf, collect(μs))
gui()

## Try eigenvalues of SSBM blocks
m = 50
k = 10
p = 0.20
#μs = linspace(0.5,0.5,1)
μs = [0.3, 0.5, 0.6]
for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))
  M = q*ones(k,k)
  M += (p-q)*eye(k)
  #M1 = [p (k-1)*q; (k-1)*q (k*p + (k-1)*(k-2)*q)]
  #d1 = vec(sum(M1,2))
  #W1 = (M1./sqrt(d))./sqrt(d')
  #@show 1.0-eigvals(W1)
  d = vec(sum(M,2))
  L = speye(k) - (M./sqrt(d))./sqrt(d')

  #@show eigvals(L)
  @show μ, ((k-1)*q/(p+(k-1)*q))
  #@show eigvals((M*M).*M)
end

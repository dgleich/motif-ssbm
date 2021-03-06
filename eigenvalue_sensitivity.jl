# In this one, we are looking at various conditioning and sensitivity measures
# of the problem.

## Setup
include("motif_codes.jl")
using MatrixNetworks
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false" # turn off Atom plotpane
using Plots
using Loess
using LaTeXStrings

## Conditioning measures and a sweep across $\mu$

# Two little helpers
zero_thresh(x, t) = x[x .>= t]
abszero_thresh(x, t) = x[abs(x) .>= t]

"""
return a list of disconnected pieces
"""
function disconnected_nontrivial_pieces(A, x)
  T = typeof(A)
  map = strong_components_map(A)
  rval = Vector{T}()
  rval2 = Vector{typeof(x)}()
  for c = 1:maximum(map)
    filt = map .== c
    if sum(filt) > 1
      push!(rval, A[filt,filt])
      push!(rval2, x[filt])
    end
  end
  return rval, rval2
end

disconnected_nontrivial_pieces(A) = disconnected_nontrivial_pieces(A,zeros(size(A,1)))[1]

function pseudo_cond(A)
  vals = svdvals(A)
  vals = vals[vals .>= 10*eps(1.0)]
  return maximum(vals)/minimum(vals)
end


function disconnected_condition_num(A)
  return maximum(map(X -> pseudo_cond(full(X)), disconnected_nontrivial_pieces(A)))
end

function connected_fiedler_condition_num(Acc)
  f,lam2 = fiedler_vector(Acc)
  d = vec(sum(Acc,1))
  z = d.*f
  return norm(f)*norm(z)/abs(dot(f,z)), lam2
end

function disconnected_fiedler_condition_num(A)
  rvals = map(X -> connected_fiedler_condition_num(X), disconnected_nontrivial_pieces(A))
  return (maximum(map(x -> x[1], rvals)),
          minimum(map(x -> x[2], rvals)))
end

function nonzero_eigenvector_condition_num(A, d::Any=1.0)
  vals,V = eig(full(A))
  X = V[:,find(abs(vals) .>= 10*eps(1.0))] # get the non-null eigenvectors
  return cond(d.*X) # row-wise scaling due to column-repetition
end


""" Get the gap between the k and k+1st non-zero normalized Laplacian eignenvalues. """
function fiedler_gap(A,k; ratio=false)
  if ratio
    gapfun = x -> begin
      if k==1
        return (2-x[1])/2.0
      else
        return (2-x[k])/(2-x[k-1])
      end
    end
  else
    gapfun = x -> begin
      if k==1
        return x[1]
      else
        return x[k]-x[k-1]
      end
    end
  end
  mats = Motif.network_matrices(A)[1]
  L = mats[:NormalizedLaplacian]
  vals = Float64[]
  gaps = Float64[]
  for Lcc = disconnected_nontrivial_pieces(L)
    Lcc_vals = eigvals!(full(L))[2:end]
    push!(gaps, gapfun(Lcc_vals))
    push!(vals, Lcc_vals...)
  end
  sort!(vals)
  return gapfun(vals), minimum(gaps), vals
end

function conditioning_measures(mat)
  mats,ops = Motif.network_matrices(mat)
  dhalf = sqrt(vec(sum(mat,1)))

  vals = Dict{Symbol,Float64}()

  vals[:NormalizedLaplacian] = disconnected_condition_num(mats[:NormalizedLaplacian])
  #vals[:NormalizedLaplacian] = cond(full(mats[:NormalizedLaplacian]))
  vals[:CombinatorialLaplacian] = disconnected_condition_num(mats[:CombinatorialLaplacian])
  vals[:WeightedAdjacency] = cond(full(mats[:WeightedAdjacency]))
  vals[:Adjacency] = cond(full(mats[:WeightedAdjacency]))

  kappaf, lam2 = disconnected_fiedler_condition_num(mat)
  vals[:Fiedler] = kappaf
  vals[:Lambda2] = lam2
#
#  vals[:NormalizedLaplacianEigenspace] = maximum(
#    map((x,idhalf) -> nonzero_eigenvector_condition_num(x,idhalf),
#    disconnected_nontrivial_pieces(mats[:NormalizedLaplacian], 1./dhalf)))

  gap,cgap,evals = fiedler_gap(mat, 2; ratio=true)
  vals[:ShiftedLaplacianRatio] = gap
  vals[:InternalShiftedLaplacianRatio] = cgap

  gap,cgap = fiedler_gap(mat, 2)
  vals[:ShiftedLaplacianGap] = gap
  vals[:InternalShiftedLaplacianGap] = cgap

  vals[:NormalizedLambda2] = evals[1]
  vals[:NormalizedLambda3] = evals[2]
  vals[:NormalizedLambda4] = evals[3]
  vals[:NormalizedLambda5] = evals[4]
  vals[:NormalizedLambda6] = evals[5]
  vals[:NormalizedLambda7] = evals[6]
  vals[:NormalizedLambda8] = evals[7]
  vals[:NormalizedLambda9] = evals[8]
  vals[:NormalizedLambda10] = evals[9]
  vals[:NormalizedLambda11] = evals[10]
  vals[:NormalizedLambda12] = evals[11]
  vals[:NormalizedLambda13] = evals[12]

  vals[:ShiftedLaplacianRatio10] = evals[9]/evals[10]

  return vals
end

m = 50
k = 10
p = 0.2
q = 0.01
A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
@show conditioning_measures(A)
##
m = 50
k = 10
p = 0.2
srand(1)
npts = 100
μs = linspace(0.1,0.8,npts)
#ntrials = 1  # just use more points instead
muvals = []

vals_A = []
vals_M = []

for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))

  A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
  M = (A*A).*A

  push!(vals_A, conditioning_measures(A))
  push!(vals_M, conditioning_measures(M))
end


##
pyplot(size=(400,300))
h = scatter(μs, map(x -> x[:ShiftedLaplacianRatio], vals_A),  label="Edge", xlim=(0.1,0.8))
scatter!(μs, map(x -> x[:ShiftedLaplacianRatio], vals_M), label="Triangle")
# we need Loess smoothed estimates of these
model = loess(collect(μs), collect(Float64, map(x -> x[:ShiftedLaplacianRatio], vals_A)))
plot!(μs, predict(model, μs), label="Edge", color=1, label="")
model = loess(collect(μs), collect(Float64, map(x -> x[:ShiftedLaplacianRatio], vals_M)))
plot!(μs, predict(model, μs), label="Triangle", color=2, label="", grid=false)
xlabel!( "μ"  )
ylabel!(L"$λ_3/λ_2$")
gui()
savefig("eigenvalue-ratio.pdf")

##


##
pyplot(size=(400,300))
h = scatter(μs, map(x -> x[:ShiftedLaplacianGap], vals_A),  label="Edge", xlim=(0.1,0.8))
scatter!(μs, map(x -> x[:ShiftedLaplacianGap], vals_M), label="Triangle")
# we need Loess smoothed estimates of these
model = loess(collect(μs), collect(Float64, map(x -> x[:ShiftedLaplacianGap], vals_A)))
plot!(μs, predict(model, μs), label="Edge", color=1, label="")
model = loess(collect(μs), collect(Float64, map(x -> x[:ShiftedLaplacianGap], vals_M)))
plot!(μs, predict(model, μs), label="Triangle", color=2, label="", grid=false)
xlabel!( "μ"  )
ylabel!(L"$λ_2-λ_3$")
gui()
savefig("eigenvalue-gap.pdf")


##
pyplot(size=(400,300))
h = scatter(μs, map(x -> x[:ShiftedLaplacianRatio10], vals_A),  label="Edge", xlim=(0.1,0.8))
scatter!(μs, map(x -> x[:ShiftedLaplacianRatio10], vals_M), label="Triangle", grid=false)
xlabel!( "μ"  )
ylabel!(L"$λ_{10}/λ_{11}$")
gui()
savefig("eigenvalue-ratio10.pdf")

##
pyplot(size=(400,300))
h = scatter(μs, map(x -> x[:ShiftedLaplacianRatio10], vals_A),  label="Edge", xlim=(0.1,0.8))
scatter!(μs, map(x -> x[:ShiftedLaplacianRatio10], vals_M), label="Triangle", grid=false)
xlabel!( "μ"  )
ylabel!(L"$λ_{10}/λ_{11}$")
gui()
savefig("eigenvalue-ratio10.pdf")

##
scatter(μs, map(x -> x[:Lambda2], vals_A),  label="Edge", xlim=(0.1,0.8))
scatter!(μs, map(x -> x[:Lambda2], vals_M), label="Triangle", grid=false)
xlabel!( "μ"  )
ylabel!(L"$λ_{2}$")
gui()
savefig("eigenvalue-lambda2.pdf")

##
scatter(μs, map(x -> x[:Fiedler], vals_A),  label="Edge", xlim=(0.1,0.8))
scatter!(μs, map(x -> x[:Fiedler], vals_M), label="Triangle", grid=false)
xlabel!( "μ"  )
ylabel!("Eigenvector conditioning\n"*L"$1/\cos(x,Dx)$")
gui()
savefig("fiedler-conditioning.pdf")


## We didn't end up using these plots

##
scatter(μs, map(x -> x[:NormalizedLaplacian], vals_A), yscale=:log10, label="Edge")
scatter!(μs, map(x -> x[:NormalizedLaplacian], vals_M), label="Triangle")
gui()
##
scatter(μs, map(x -> x[:Lambda2], vals_A),  label="Edge")
scatter!(μs, map(x -> x[:Lambda2], vals_M), label="Triangle")
gui()
##
scatter(μs, map(x -> x[:Fiedler], vals_A),  label="Edge")
scatter!(μs, map(x -> x[:Fiedler], vals_M), label="Triangle")
gui()

##
scatter(μs, map(x -> x[:NormalizedSimpleSep], vals_A),  label="Edge", yscale=:log10)
scatter!(μs, map(x -> x[:NormalizedSimpleSep], vals_M), label="Triangle")
gui()


##
# Look at the smallest few eigenvalues
plotkeys = [:NormalizedLambda2, :NormalizedLambda3, :NormalizedLambda4, :NormalizedLambda5,
            :NormalizedLambda6, :NormalizedLambda7, :NormalizedLambda8, :NormalizedLambda9,
            :NormalizedLambda10, :NormalizedLambda11,  ]
plot() # initialize a plot
for k in plotkeys
  vals = map(x -> x[k], vals_M)
  scatter!(μs, vals)
end
gui()

##
plot() # initialize a plot
for k in plotkeys
  vals = map(x -> x[k], vals_A)
  scatter!(μs, vals)
end
gui()

##
# Okay, now look at the ratio at 10
scatter(μs, map(x -> x[:ShiftedLaplacianRatio10], vals_A),  label="Edge")
scatter!(μs, map(x -> x[:ShiftedLaplacianRatio10], vals_M), label="Triangle")
# we need Loess smoothed estimates of these
#model = loess(collect(μs), collect(Float64, map(x -> x[:ShiftedLaplacianRatio10], vals_A)))
#plot!(μs, predict(model, μs), label="Edge")
#model = loess(collect(μs), collect(Float64, map(x -> x[:ShiftedLaplacianRatio10], vals_M)))
#plot!(μs, predict(model, μs), label="Triangle")
gui()

## One of the questions we had to resolve was whether the gap lambda2
# is smaller for motifs vs. matrices. FOr this, we need a slightly more refined analysis.
m = 50
k = 10
p = 0.2
srand(1)
npts = 10
μs = linspace(0.1,0.8,npts)
ntrials = 20

gap_larger = zeros(npts, ntrials)
for (i,μ) in enumerate(μs)
  for t = 1:ntrials
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A

    gapA = fiedler_gap(A,2; ratio=true)[1]
    gapM = fiedler_gap(M,2; ratio=true)[1]

    if gapA >= gapM
      gap_larger[i,t] = 1.0
    end
  end
end
plot(μs, mean(gap_larger,2))
gui()

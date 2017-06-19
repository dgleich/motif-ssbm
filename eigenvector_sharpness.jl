# In this one, we are looking at various conditioning and sensitivity measures
# of the problem.

## Setup
include("motif_codes.jl")
using MatrixNetworks
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false" # turn off Atom plotpane
using Plots
using Loess
using LaTeXStrings


## Let's just show a random combination of eigenvectors
m = 50
k = 10
p = 0.2
srand(1)
μ = 0.41
A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
M = (A*A).*A
q = μ*p/((1-μ)*(k-1))

LA = Motif.network_matrices(A)[1][:NormalizedLaplacian]
LM = Motif.network_matrices(M)[1][:NormalizedLaplacian]

XA = eigvecs(collect(LA))[:,2:10]
XM = eigvecs(collect(LM))[:,2:10]

pyplot(size=(1200,600))
C = randn(9,1) # 3 random combinations
C = C ./ sqrt(sum(C.^2, 1))
plot(scatter(XA*C), scatter(XM*C), leg=false)
gui()



##
# Show the eigenvectors
plot(scatter(XA[:,1:3]), scatter(XM[:,1:3]), leg=false)
gui()
##
# Show the eigenctors with smoothing
pyplot(size=(400,400))
plot(XA[:,1:3], xticks=0:50:500,yticks=[0],legend=false,background=false,border=:black,linewidth=0.5)
for i=1:3
  model = loess(collect(1.0:m*k), XA[:,i], span=0.1)
  plot!(predict(model, collect(1.0:m*k)), color=i)
end
gui()
savefig("eigenvector-localization-edges-$(m)-$(k)-$(p)-$(μ).pdf")

##
plot(XM[:,1:3], xticks=0:50:500,yticks=[0],legend=false,background=false,border=:black,linewidth=0.5)
for i=1:3
  model = loess(collect(1.0:m*k), XM[:,i], span=0.1)
  plot!(predict(model, collect(1.0:m*k)), color=i)
end
gui()
savefig("eigenvector-localization-motif-$(m)-$(k)-$(p)-$(μ).pdf")


## This shows mostly what we want.

## Let's do the test on random recovery
m = 50
k = 10
p = 0.2
srand(1)

M = (A*A).*A

npts = 100
μs = linspace(0.1,0.8,npts)

ntrials = 50

RA = zeros(npts, ntrials)
RM = zeros(npts, ntrials)

scorevec = x -> begin
  perm = sortperm(vec(x))
  return max(Motif.score_set(perm[1:m],m,k), Motif.score_set(perm[end-m+1:end],m,k))
end

for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))

  A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
  M = (A*A).*A

  dahalf = vec(sqrt(sum(A,1)))
  dmhalf = vec(sqrt(sum(M,1)))

  LA = Motif.network_matrices(A)[1][:NormalizedLaplacian]
  LM = Motif.network_matrices(M)[1][:NormalizedLaplacian]

  XA = eigvecs(collect(LA))[:,2:10]
  XM = eigvecs(collect(LM))[:,2:10]

  for t=1:ntrials

    C = randn(9,1) # 3 random combinations
    C = C ./ sqrt(sum(C.^2, 1))

    xa = XA*C
    xm = XM*C

    RA[i,t] = scorevec(xa./dahalf)
    RM[i,t] = scorevec(xm./dahalf)

  end
end

##
plot(μs, median(RA,2))
plot!(μs, median(RM,2))
gui()


## These were some initial localizatin studies.

m = 50
k = 10
p = 0.2
srand(1)
npts = 20
μs = linspace(0.2,0.4,npts)
#ntrials = 1  # just use more points instead
muvals = []

vals_A = []
vals_M = []

locA = []
locM = []

pyplot(size=(1600,800))

for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))

  A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
  M = (A*A).*A

  LA = Motif.network_matrices(A)[1][:NormalizedLaplacian]
  LM = Motif.network_matrices(M)[1][:NormalizedLaplacian]

  XA = eigvecs(collect(LA))[:,2:10]
  XM = eigvecs(collect(LM))[:,2:10]

  prA = 1.0./sum(XA.^4,1)
  prM = 1.0./sum(XM.^4,1)

  push!(locA, prA)
  push!(locM, prM)


end

plot(
  heatmap(vcat(locA...), legend=false, yticks=[]),
  heatmap(vcat(locM...), legend=false, yticks=[])
)
gui()
##

plot(μs, sum(vcat(locA...),2))
plot!(μs, sum(vcat(locM...),2))
gui()
##

for (i,μ) in enumerate(μs)
  q = μ*p/((1-μ)*(k-1))

  A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
  M = (A*A).*A

  LA = Motif.network_matrices(A)[1][:NormalizedLaplacian]
  LM = Motif.network_matrices(M)[1][:NormalizedLaplacian]

  XA = eigvecs(collect(LA))[:,2:4]
  XM = eigvecs(collect(LM))[:,2:4]

  prA = 1.0./sum(XA.^4,1)
  prM = 1.0./sum(XM.^4,1)
  @show prA, sum(prA)
  @show prM, sum(prM)


  plot(
    plot(XA, legend=false, yticks=[], ylim=[-0.3,0.3]),
    plot(XM, legend=false, yticks=[], ylim=[-0.3,0.3])
  )
  gui()
  sleep(0.5)

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

## Setup
include("motif_codes.jl")
using MatrixNetworks
using Plots


## Power method accuracy
# We are going to look at accuracy in the power method as a function of
# * number of steps
# * mu
# * ntrials


# Here is one experiment for a given matrix
function step_vecs_experiment(mat, m, k, nsteps)
  mats,ops = Motif.network_matrices(mat)
  dhalf = sqrt(vec(sum(mat,1)))
  S = ops[:ShiftedLaplacian]

  x0 = ones(size(mat,1))

  nmax = maximum(nsteps)
  result = zeros(length(nsteps))
  xstep = x0
  for step=1:nmax
    xstep = S*xstep
    xstep = xstep/norm(xstep)

    for j in findin(nsteps, step)
      # hopefully just one, so this shouldn't duplicate work, but we keep it simple!
      perm = sortperm(vec(xstep)./dhalf)
      result[j] = max(Motif.score_set(perm[1:m],m,k), Motif.score_set(perm[end-m+1:end],m,k))
    end
  end
  return result
end

srand(2)
m = 50
k = 10
p = 0.20

npts = 50
μs = linspace(0.1,0.8,npts)
ntrials = 200
nsteps = vec([1 2 3 4 5 7 10 15 20 30 40 50 75 100 150 200 300 400])

accs_edge = zeros(npts, ntrials, length(nsteps))
accs_motif = zeros(npts, ntrials, length(nsteps))


# Here is the experiment sweep
for (i,μ) = enumerate(μs)
  q = μ*p/((1-μ)*(k-1))
  for t = 1:ntrials
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A+A

    accs_edge[i,t,:] = step_vecs_experiment(A, m, k, nsteps)
    accs_motif[i,t,:] = step_vecs_experiment(M, m, k, nsteps)
  end
end

## Final plot
pyplot(size=(500,300))
l = @layout([a b]) # give a more width to account for colorbar
plt = plot(layout=l)
contour!(plt[1], map(x -> string(x), nsteps), # map matvecs to string
  collect(μs), squeeze(mean(accs_edge, 2),2),
  xlabel="matrix-vector products",
  ylabel="μ",
  clims=(0.15,1), fill=true, levels=5, colorbar=false)
contour!(plt[2], map(x -> string(x), nsteps),
  collect(μs), squeeze(mean(accs_motif, 2),2),
  xlabel="matrix-vector products",
  ylabel=" ",
  clims=(0.15,1), fill=true, levels=5)
xaxis!(rotation=45)
gui()
savefig("powermethod-accuracy.pdf")

## Make a super-high matvec plot
# What happens in the previous picture if we go for thousands of iterations?
# Turns out -- nothing! The plot stays the same, it just gets longer.
# An interesting case is mu=0.6, where we still get good accuracy for motifs
# but bad accuracy for edges.

npts = 3
μs = linspace(0.55,0.65,npts)
ntrials = 200
nsteps = vec([1 2 3 4 5 7 10 15 20 30 40 50 75 100 150 200 300 400 500 700 1000 1500 2000 2500])

accs_edge = zeros(npts, ntrials, length(nsteps))
accs_motif = zeros(npts, ntrials, length(nsteps))


# Here is the experiment sweep
for (i,μ) = enumerate(μs)
  q = μ*p/((1-μ)*(k-1))
  for t = 1:ntrials
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A+A

    accs_edge[i,t,:] = step_vecs_experiment(A, m, k, nsteps)
    accs_motif[i,t,:] = step_vecs_experiment(M, m, k, nsteps)
  end
end

##
l = @layout([a b]) # give a more width to account for colorbar
plt = plot(layout=l)
contour!(plt[1], map(x -> string(x), nsteps), # map matvecs to string
  collect(μs), squeeze(mean(accs_edge, 2),2),
  xlabel="matrix-vector products",
  ylabel="μ",
  clims=(0.15,1), fill=true, levels=5, colorbar=false)
contour!(plt[2], map(x -> string(x), nsteps),
  collect(μs), squeeze(mean(accs_motif, 2),2),
  xlabel="matrix-vector products",
  ylabel=" ",
  clims=(0.15,1), fill=true, levels=5)
#


## Other possible plots
l = @layout([a b; c d])
plt = plot(layout=l, legend=false)
heatmap!(plt[1], map(x -> string(x), nsteps), collect(μs), squeeze(mean(accs_edge, 2),2), clims=(0.15,1))
heatmap!(plt[3], map(x -> string(x), nsteps), μs, squeeze(mean(accs_motif, 2),2), clims=(0.15,1))
contour!(plt[2], squeeze(mean(accs_edge, 2),2), clims=(0.15,1), fill=true, levels=5)
contour!(plt[4], squeeze(mean(accs_motif, 2),2), clims=(0.15,1), fill=true, levels=5)
gui()
#@show mean(accs_motif, 2)
##
contour(squeeze(mean(accs_edge, 2),2), clims=(0.15,1), fill=true, levels=5)
gui()

##
contour(squeeze(mean(accs_motif, 2),2), clims=(0.15,1), fill=true, levels=5)
gui()

##
l = @layout([a b])
plt = plot(layout=l, legend=false)
contour!(plt[1], map(x -> string(x), nsteps), collect(μs), squeeze(mean(accs_edge, 2),2), clims=(0.15,1), fill=true, levels=5)
contour!(plt[2],  map(x -> string(x), nsteps), collect(μs), squeeze(mean(accs_motif, 2),2), clims=(0.15,1), fill=true, levels=5)
gui()

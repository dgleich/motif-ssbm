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
μs = linspace(0.4,0.8,npts)
ntrials = 100
nsteps = vec([1 2 3 4 5 7 10 15 20 30 40 50 75 100 150 200 300 400])

accs_tri = zeros(npts, ntrials, length(nsteps))
accs_clique = zeros(npts, ntrials, length(nsteps))


# Here is the experiment sweep
for (i,μ) = enumerate(μs)
  q = μ*p/((1-μ)*(k-1))
  for t = 1:ntrials
    A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
    M = (A*A).*A
    W = Motif.clique4_weighted(A)+M

    accs_tri[i,t,:] = step_vecs_experiment(M, m, k, nsteps)
    accs_clique[i,t,:] = step_vecs_experiment(W, m, k, nsteps)
  end
end

## Final plot
pyplot(size=(500,300))
l = @layout([a b]) # give a more width to account for colorbar
plt = plot(layout=l)
contour!(plt[1], map(x -> string(x), nsteps), # map matvecs to string
  collect(μs), squeeze(mean(accs_tri, 2),2),
  xlabel="matrix-vector products",
  ylabel="μ",
  clims=(0.15,1), fill=true, levels=5, colorbar=false)
contour!(plt[2], map(x -> string(x), nsteps),
  collect(μs), squeeze(mean(accs_clique, 2),2),
  xlabel="matrix-vector products",
  ylabel=" ",
  clims=(0.15,1), fill=true, levels=5)
xaxis!(rotation=45)
hline!(plt[1], [detect_thresh_mu(m,k,p)], legend=false, color=:lightgreen)
hline!(plt[2], [detect_thresh_mu(m,k,p)], legend=false, color=:lightgreen)
gui()
savefig("powermethod-accuracy-clique4.pdf")

## Setup
include("motif_codes.jl")
using MatrixNetworks
using Plots


## Power method animation
# In power method steps, we saw that μ=0.6 gave good results for enough matvecs for
# the motif matrix, but bad results for the edge case. (or not as clean).
# Here, we'd like to animate these matvecs to see what's going on.

function animate_powermethod(mat, nmax, gifname )
    mats,ops = Motif.network_matrices(mat)
    dhalf = sqrt(vec(sum(mat,1)))
    S = ops[:ShiftedLaplacian]

    x0 = ones(size(mat,1))

    xstep = x0

    anim = @animate for step=1:nmax
      xstep = S*xstep
      xstep = xstep/norm(xstep)

      l = @layout([a;b])
      plt = plot(layout=l)
      z = xstep./dhalf
      scatter!(plt[1], z, ylim=(-0.035,0.035), ylabel="Iteration $(step)",
        legend=false, border=false, left_margin=10px)

      perm = sortperm(z)
      bottomset = Motif.best_set_and_score(perm[1:m],m,k)
      topset = Motif.best_set_and_score(perm[end-m+1:end],m,k)


      flipsign = 1
      if topset[2] >= bottomset[2]
        bset = topset[1]
      else
        bset = bottomset[1]
        flipsign = -1
      end


      z1 = minimum(flipsign*z[(m*(bset-1)+1):(m*bset)])
      @show z1

      hline!(plt[1], [flipsign*z1])
      correct = intersect(find(flipsign*z .>= z1), (m*(bset-1)+1):(m*bset))
      errs = find(flipsign*z .>= z1)
      errs = setdiff(errs, (m*(bset-1)+1):(m*bset)) # no errors in 1:m
      scatter!(plt[1], errs, z[errs], color=2)
      scatter!(plt[1], correct, z[correct], color=3)

      scatter!(plt[2], sort(xstep./dhalf), ylim=(-0.035,0.035),
        legend=false, border=false )
      xticks!(0:50:500)
      yticks!(Int[])
      gui()
      sleep(0.01)
  end
  gif(anim, gifname, fps=20)
end

m = 50
k = 10
p = 0.20
μ = 0.41
maxstep=200
q = μ*p/((1-μ)*(k-1))

srand(1)
A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
M = (A*A).*(A)
pyplot(size=(400,600))

animate_powermethod(M, maxstep, "power-motif-mu=$(μ).gif")

##
animate_powermethod(A, maxstep, "power-edges-mu=$(μ).gif")

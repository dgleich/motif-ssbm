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
      scatter!(plt[1], xstep./dhalf, ylim=(-0.035,0.035), ylabel="Iteration $(step)",
        legend=false, border=false, left_margin=10px)
      scatter!(plt[2], sort(xstep./dhalf), ylim=(-0.035,0.035),
        legend=false, border=false )
      xticks!(0:50:500)
      yticks!(Int[])
      #gui()
      #sleep(0.05)
  end
  gif(anim, gifname, fps=15)
end

m = 50
k = 10
p = 0.20
μ=0.6
maxstep=150
q = μ*p/((1-μ)*(k-1))

srand(1)
A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
M = (A*A).*(A)
pyplot(size=(400,600))
animate_powermethod(A, maxstep, "edges-mu=0.6.gif")

##
animate_powermethod(M, maxstep, "motif-mu=0.6.gif")

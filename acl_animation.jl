## Setup
include("motif_codes.jl")
using MatrixNetworks
using DataStructures
using Plots


## ACL method animation
# In power method steps, we saw that μ=0.6 gave good results for enough matvecs for
# the motif matrix, but bad results for the edge case. (or not as clean).
# Here, we'd like to animate these matvecs to see what's going on.

function animate_acl(mat, gifname, ylim=(0.000001,0.001))

  pattern = [100 125 150 175 200 250 300 400 500 600 700 800 900]
  nsteps = [10 15 20 30 40 50 60 70 80 90 pattern 10*pattern]
  #nsteps = [10 15 20 30 40 50 60 70 80 90]

  alpha = 0.99
  eps = 1.0e-6
  seed = 2

  A = mat
  dvec = vec(sum(A,1)) # get the degrees

  colptr = A.colptr
  rowval = A.rowval
  nzval = A.nzval

  n = size(A,1)

  x = zeros(n)
  r = zeros(n)
  Q = Queue(Int)

  r[seed] = 1.0
  enqueue!(Q,seed)

  pushcount = 0
  pushvol = 0
  stepcheck = 1

  anim = Animation()
  while length(Q) > 0 && stepcheck <= length(nsteps)
      pushcount += 1
      u = dequeue!(Q)

      du = dvec[u] # get the degree
      idu = 1.0/du

      pushval = r[u]
      x[u] = get(x,u,0.0) + (1-alpha)*pushval
      r[u] = 0.0

      pushval = pushval*alpha

      for nzi in colptr[u]:(colptr[u+1] - 1)
          pushvol += 1
          @inbounds v = rowval[nzi]
          @inbounds dv = dvec[v] # degree of v
          val = pushval*nzval[nzi]*idu

          rvold = r[v]
          rvnew = rvold + val
          r[v] = rvnew
          if rvnew > eps*dv && rvold <= eps*dv
              enqueue!(Q,v)
          end
      end

      if pushcount == nsteps[stepcheck]


        # show errors in red
        begin # use extra variables here
          z = x./dvec
          scatter(z, ylim=ylim, ylabel="Iteration $(step)",
            legend=false, border=false, left_margin=10px, yscale=:log10)
          z1 = minimum(z[find(z[1:m])])
          hline!([z1])
          errs = find(z .>= z1)
          errs = setdiff(errs, 1:m) # no errors in 1:m
          scatter!(errs, z[errs], color="red")
        end
        #scatter!(plt[2], sort(x./dvec), ylim=(0,0.01),
        #  legend=false, border=false )
        xticks!(0:50:500)
        yticks!(Int[])
        title!("step = $pushcount")
        gui()
        sleep(0.01)
        frame(anim)
        stepcheck += 1
      end
  end


  gif(anim, gifname, fps=15)
end

m = 50
k = 10
p = 0.50
μ=0.40
q = μ*p/((1-μ)*(k-1))

srand(1)
A,sets = Motif.symmetric_stochastic_block_model(m,k,p,q)
M = (A*A).*(A)
W = Motif.clique4_weighted(A)

pyplot(size=(400,300))
animate_acl(M,  "acl-motif-mu=0.4.gif", (0.00000001,0.0001))


##
animate_acl(A, "acl-edges-mu=0.4.gif")

##
animate_acl(W, "acl-clique4-mu=0.5.gif", (0.00000001,0.00001))

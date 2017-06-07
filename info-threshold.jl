##
using Roots

## Find the info theory threshold
# This is based on Abbe's SSBM paper, remark 4 on the version dated march 21, 2017.
# sbm_jmlr_4.pdf
# This says that.
function thresh(m,k,p,q)
  n = m*k
  a = p*(n)/log(n)
  b = q*(n)/log(n)
  t = (sqrt(a)-sqrt(b))^2/k
end

thresh_mu(m,k,p,μ) = thresh(m,k,p,μ*p/((1-μ)*(k-1)))

# Do bisection search to find the threshold in terms of our mixing parameter mu
function thresh_mu(m,k,p)
  fzero(mu -> thresh_mu(m,k,p,mu)-1, [0.0, (k-1)/k])
end

thresh_mu(m,k,0.5)


##
thresh_mu(m,k,0.2)

##
# The above results were for full recovery
# There is another threshold for detection
function detect_thresh_snr(m,k,p,q)
  n = m*k
  a = p*n
  b = q*n
  snr = (a-b)^2/(k*(a + (k-1)*b))
end

detect_thresh_snr_mu(m,k,p,μ) = detect_thresh_snr(m,k,p,μ*p/((1-μ)*(k-1)))

function detect_thresh_mu(m,k,p)
  fzero(mu -> detect_thresh_snr_mu(m,k,p,mu)-1, [0.0, (k-1)/k])
end
@show detect_thresh_mu(50,10,0.5)
@show detect_thresh_mu(50,10,0.2)

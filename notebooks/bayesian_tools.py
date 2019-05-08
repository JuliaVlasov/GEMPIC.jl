# ---
# jupyter:
#   jupytext:
#     comment_magics: false
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.1.1
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# %%
using BayesianTools.ProductDistributions
p = ProductDistribution(Normal(0,1), Beta(1.,1.))
n = length(p) ## 2 -> Number of distributions in the product

# %%
using Random
rand!(p, zeros(Float64,(2,100)))

# %%
using BayesianTools.Links
function mcmc_wrong(iters)
   chain = zeros(Float64, iters)
   gamma = Gamma(2, 1)
   d = Improper(0, +Inf)
   lx  = 1.0
   for i in 1:iters
      xs = link(d, lx) + randn()
      lxs = invlink(d, xs)
      a = logpdf(gamma, lxs)-logpdf(gamma, lx)       
      (rand() < exp(a)) && (lx = lxs)
      chain[i] = lx
   end
   return chain
end

# %%
function mcmc_right(iters)
   chain = zeros(Float64,iters)
   gamma = Gamma(2, 1)
   d = Improper(0, +Inf)
   lx  = 1.0
   for i in 1:iters
      xs = link(d, lx) + randn()
      lxs = invlink(d, xs)
      a = logpdf(gamma, lxs)-logpdf(gamma, lx)
      ## Log absolute jacobian adjustment
      a = a - logjacobian(d, lxs) + logjacobian(d, lx)
      (rand() < exp(a)) && (lx = lxs)
      chain[i] = lx
   end
   return chain
end

# %%
Plots.histogram([mc0, mc1], normalize=true, bins = 100, fill=:slategray, layout = (1,2), lab = "draws")



# %%

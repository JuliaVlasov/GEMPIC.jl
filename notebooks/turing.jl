# ---
# jupyter:
#   jupytext:
#     comment_magics: false
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.1
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

using StatsPlots

using Turing

# Define a simple Normal model with unknown mean and variance.
@model gdemo(x, y) = begin
  x ~ Normal()
  y ~ Normal()
  return x, y
end

#  Run sampler, collect results
chn = sample(gdemo(), SMC(10000))

# Summarise results (currently requires the master branch from MCMCChains)
describe(chn)

# Plot and save results
plot(chn)



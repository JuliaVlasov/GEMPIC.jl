# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     comment_magics: false
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Julia 1.4.0
#     language: julia
#     name: julia-1.4
# ---

# ## Initialize a distribution

using GEMPIC

@doc SumCosGaussian

# +
using Plots

k, α = 0.5, 0.1
nx, nv = 64, 128
xmin, xmax = 0, 2π/k
vmin, vmax = -6, 6
xg = LinRange(xmin, xmax, nx+1)[1:end-1] 
vg = LinRange(vmin, vmax, nv+1)[1:end-1];

# -

df_1d1v = SumCosGaussian{1,1}( [[k]] , [α] , [[1.0]] , [[0.0]], [1.0])
f = zeros(Float64, (nx,nv))
for  j in eachindex(vg), i in eachindex(xg)
    f[i,j] = df_1d1v( xg[i], vg[j])
end
surface(xg, vg, f')

two_gaussians = SumCosGaussian{1,1}( [[k]] , [0.0] , [[0.2],[0.2]] , [[-1.0],[1.0]], [0.5, 0.5])
f = zeros(Float64, (nx,nv))
for  j in eachindex(vg), i in eachindex(xg)
    f[i,j] = two_gaussians( xg[i], vg[j])
end
surface(xg, vg, f')

# +
fxv(x, v) = ( 1 + α * cos(k * x )) / sqrt(2π) * exp(- (v^2)/ 2)

surface(xg, vg, fxv)
# -

@doc CosSumGaussian

df = CosSumGaussian{1,1}( [[k]], [0.1], [[1.0]], [[0.0]], [1.0] )

plot( xg, [GEMPIC.eval_x_density(df, x) for x in xg])

plot( vg, [GEMPIC.eval_v_density(df, v) for v in vg])

# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     comment_magics: false
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.1.1
#     language: julia
#     name: julia-1.1
# ---

# ## Initialize a probablity distribution

include("../src/distributions.jl")

?SumCosGaussian

# +
using Plots

k, α = 0.5, 0.1
nx, nv = 64, 128
xmin, xmax = 0, 2π/k
vmin, vmax = -6, 6
xg = range(xmin, stop=xmax, length=nx+1)[1:end-1] |> collect
vg = range(vmin, stop=vmax, length=nv+1)[1:end-1] |> collect;

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

?CosSumGaussian

df = CosSumGaussian{1,1}( [[k]], [0.1], [[1.0]], [[0.0]], [1.0] )

plot( xg, [eval_x_density(df, x) for x in xg])

plot( vg, [eval_v_density(df, v) for v in vg])



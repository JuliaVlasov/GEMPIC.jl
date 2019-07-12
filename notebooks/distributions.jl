# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     comment_magics: false
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# +
using Plots, GEMPIC

kx, ϵ = 0.5, 0.1
nx, nv = 64, 128
xmin, xmax = 0, 2π/kx
vmin, vmax = -6, 6
xg = range(xmin, stop=xmax, length=nx+1)[1:end-1] |> collect
vg = range(vmin, stop=vmax, length=nx+1)[1:end-1] |> collect;

# +


fxv(x, v) = ( 1 + ϵ * cos(kx * x )) / sqrt(2π) * exp(- (v^2)/ 2)

surface(xg, vg, fxv)
# -

v_thermal = hcat([1.0])
v_mean    = hcat([0.0])
δ = 1.0
df = CosGaussian( (1,1), 1, 1, hcat([kx]), [0.1], v_thermal, v_mean, δ )

plot( xg, [eval_x_density(df, x) for x in xg])

plot( vg, [eval_v_density(df, v) for v in vg])

using Distributions
using Random

# $$
#  f(x,v) = \frac{1}{\sqrt{2π}} \big( 1 + ϵ \cos ( k_x x) \big) e^{-v^2/2} 
# $$

?Cosine

# $$
# f(x;\mu,\sigma)=\frac{1}{2\sigma}
# \left[1+\cos\left(\frac{x-\mu}{\sigma}\,\pi\right)\right]\,
# $$

?Normal

fx = Cosine(1,1)
fv = Normal(0,1)

plot(vg, pdf.(fv, vg))

x = rand!(fx, zeros(10000))
v = rand!(fv, zeros(10000));

Plots.histogram(x, normalize=true, bins = 100, fill=:slategray, lab = "draws")
plot!(xg, pdf.(fx, xg))



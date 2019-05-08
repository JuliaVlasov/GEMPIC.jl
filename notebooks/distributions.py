# -*- coding: utf-8 -*-
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
using Distributions
using Random

# %% [markdown]
# $$
#  f(x,v) = \frac{1}{\sqrt{2π}}( 1 + ϵ \cos ( k_x x)  e^{-v^2/2} 
# $$

# %%
?Cosine

# %% [markdown]
# $$
# f(x;\mu,\sigma)=\frac{1}{2\sigma}
# \left[1+\cos\left(\frac{x-\mu}{\sigma}\,\pi\right)\right]\,
# $$

# %%
?Normal

# %%
kₓ, ϵ = 0.5, 0.1

# %%
using Plots

# %%
nx, nv = 64, 128
xmin, xmax = 0, 2π
vmin, vmax = -6, 6
xg = range( 0, stop=2, length=nx+1)[1:end-1] |> collect
vg = range(-6, stop=6, length=nx+1)[1:end-1] |> collect;

# %%
fx = Cosine(1,1)
fv = Normal(0,1)

# %%

# %%
plot(vg, pdf.(fv, vg))

# %%
rand(fx)

# %%
x = rand!(fx, zeros(10000))
v = rand!(fv, zeros(10000));

Plots.histogram(x, normalize=true, bins = 100, fill=:slategray, lab = "draws")
plot!(xg, pdf.(fx, xg))

# %%
fxv(x, v) = ( 1 + ϵ * cos(kₓ * x )) / sqrt(2π) * exp(- (v^2)/ 2)

# %%
surface(xg, vg, fxv)

# %%

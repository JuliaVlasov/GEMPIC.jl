# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Julia 1.4.1
#     language: julia
#     name: julia-1.4
# ---

using Distributed
using BenchmarkTools
using Random
using Test
rmprocs(workers())
addprocs(4)
nworkers()

@everywhere using SharedArrays

# +
# shared array
function compute_rho( steps)
    nx = 64
    np = 1000000
    xmin, xmax = -6, 6
    Lx = xmax - xmin
    rng = MersenneTwister(123)
    xp = randn(rng, np);
    wp = similar(xp)
    fill!(wp, 1.0 / np)
    rho = SharedArray(zeros(Float64, nx))
    
    for step = 1:steps

        fill!(rho, 0.0)
        @sync @distributed (+) for i in 1:length(rho)
            x_norm = (xp[i]-xmin) / Lx
            ip = trunc(Int,  x_norm * nx)+1
            rho[ip] += wp[i]
        end

    end

    rho

end

@time compute_rho(100);
# -

@test sum(rho) ≈ 1.0



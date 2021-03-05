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

import Base.Threads: @sync, @spawn, nthreads, threadid
using Random
using Test

const nx = 64
const np = 1000000

# +
function serial_deposition( steps )

    xmin, xmax = -6, 6
    Lx = xmax - xmin
    rho = zeros(Float64, nx)
    rng = MersenneTwister(123)
    xp = randn(rng, np);
    wp = similar(xp)
    fill!(wp, 1.0 / np)
    ntid = nthreads()

    for step = 1:steps

        fill!(rho, 0.0)

        for i in 1:np
            x_norm = (xp[i]-xmin) / Lx
            ip = trunc(Int,  x_norm * nx)+1
            rho[ip] += wp[i]
        end

    end

    rho

end

# +
function parallel_deposition( steps )

    xmin, xmax = -6, 6
    Lx = xmax - xmin
    rho = zeros(Float64, nx)
    rng = MersenneTwister(123)
    xp = randn(rng, np);
    wp = similar(xp)
    fill!(wp, 1.0 / np)
    ntid = nthreads()
    rho_local = [zero(rho) for _ in 1:ntid]
    chunks = Iterators.partition(1:np, np÷ntid)

    for step in 1:steps

        @sync for chunk in chunks
            @spawn begin
                tid = threadid()
                fill!(rho_local[tid], 0.0)
                for i in chunk
                    x_norm = (xp[i]-xmin) / Lx
                    ip = trunc(Int,  x_norm * nx)+1
                    rho_local[tid][ip] += wp[i]
                end
            end
        end

        rho .= reduce(+,rho_local)

     end

     rho

end
# -

@show nthreads()

serial_deposition(1) # trigger compilation
@time rho = serial_deposition(100)
@test sum(rho) ≈ 1.0

parallel_deposition(1) # trigger compilation
@time rho = parallel_deposition(100)
@test sum(rho) ≈ 1.0

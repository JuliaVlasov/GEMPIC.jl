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

import Base.Threads: @sync, @spawn, nthreads, threadid, threads
using Random
using Plots
using BenchmarkTools

@sync for i = 1:12
    @spawn println(i, " on thread ", threadid())
end

nx = 64
x = LinRange(-6, 6, nx)
rho = zeros(Float64, nx)
np = 2000000
rng = MersenneTwister(123)
xp = randn(rng, np);

function find_index( xp :: Float64, xmin, xmax, nx) :: Int
    x_norm = (xp-xmin)/(xmax-xmin)
    trunc(Int,  x_norm * nx) + 1
end

# +
function serial_deposition!(rho, xp, x)

    fill!(rho, 0.0)
    nx = size(rho)[1]
    np = size(xp)[1]
    xmin = first(x)
    xmax = last(x)
    for i in eachindex(xp)
        ip = find_index(xp[i], xmin, xmax, nx)
        rho[ip] += 1 / np
    end
    
end
# -

@time serial_deposition!(rho, xp, x)

plot(x, rho, m=:o)

function threaded_deposition!(rho, rho_local, xp, x)
    
     fill!(rho_local, 0.0)
     nx = size(rho)[1]
     np = size(xp)[1]
     xmin, xmax = first(x), last(x)
     @threads for i = 1:np
             tid = threadid()
                 ip = find_index(xp[i], xmin, xmax, nx)
                 rho_local[ip, tid] += 1 / np
     end
    
     rho .= vec(sum(rho_local, dims=2))
end

rho_local = zeros(Float64, ( nx, nthreads()))
@time threaded_deposition!(rho, rho_local, xp, x)
@show round(sum(rho), digits=6)
plot(x, rho)

V = zeros(8)
for i = 1:100000
    cc = rand(1:8)
    V[cc] += 1
end
V

# +
using Random
using Base.Threads

V = let
    mt = Tuple([MersenneTwister() for _ in 1:nthreads()])
    Vv = Tuple([zeros(3) for _ in 1:nthreads()])
    @threads for i = 1:100000
        @inbounds cc = rand(mt[threadid()], 1:3)
        @inbounds Vv[threadid()][cc] += 1
    end
    reduce(+, Vv)
end

# +
function worker(iters, rng)
    v = zeros(3)
    for i = 1:iters
        cc = rand(rng, 1:3)
        v[cc] += 1
    end
    v
end

V = let
    mt = Tuple([MersenneTwister() for _ in 1:nthreads()])
    Vv = [zeros(3) for _ in 1:nthreads()]
    jobs_per_thread = fill(div(100000, nthreads()),nthreads())
    for i in 1:100000-sum(jobs_per_thread)
        jobs_per_thread[i] += 1
    end
    @assert sum(jobs_per_thread) == 100000
    @threads for i = 1:nthreads()
        Vv[threadid()] = worker(jobs_per_thread[threadid()], mt[threadid()])
    end
    reduce(+, Vv)
end

# +
# examples/03-reduce.jl
using MPI
MPI.Init()

comm = MPI.COMM_WORLD
root = 0
psize = MPI.Comm_size(comm)
prank = MPI.Comm_rank(comm)

rng = MersenneTwister()
xp = randn(rng, 1000000)
Iterators.partition(1:np, nthreads())

r = 

sr = MPI.Reduce(r, +, root, comm)

if MPI.Comm_rank(comm) == root
    println("sum of ranks = $sr")
end

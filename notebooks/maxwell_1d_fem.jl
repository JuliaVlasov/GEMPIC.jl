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
#     display_name: Julia 1.2.0
#     language: julia
#     name: julia-1.2
# ---

# # Maxwell solver FEM 1D

using Plots

const mode = 2

include("../src/low_level_bsplines.jl")
include("../src/maxwell_1d_fem.jl")

eta1_min = .0
eta1_max = 2π
nc_eta1  = 128
Lx = eta1_max - eta1_min
delta_eta1 = Lx / nc_eta1
domain = [eta1_min, eta1_max]
deg = 3

maxwell_1d = Maxwell1DFEM(domain, nc_eta1, deg);

ex = zeros(Float64, nc_eta1)
ey = zeros(Float64, nc_eta1)
bz = zeros(Float64, nc_eta1);

bz_exact = zeros(Float64, nc_eta1)
ex_exact = zeros(Float64, nc_eta1)
ey_exact = zeros(Float64, nc_eta1)
rho      = zeros(Float64, nc_eta1)
sval     = zeros(Float64, nc_eta1);

cos_k(x) = cos(mode*2*pi*x/Lx) 
sin_k(x) = sin(mode*2*pi*x/Lx) 

# Test Poisson
#-------------
# Set exact solution
for i = 1:nc_eta1
   xi = eta1_min + (i-1)*delta_eta1
   ex_exact[i] = sin_k(xi)/(2.0*mode*pi/Lx)
end

x = range(eta1_min , stop=eta1_max, length=nc_eta1+1)[1:end-1] |> collect;

compute_rhs_from_function!( rho, maxwell_1d, cos_k, deg)

compute_e_from_rho!( ex, maxwell_1d, rho ) 

plot(x, ex_exact)
sval = eval_uniform_periodic_spline_curve(deg-1, ex);

scatter!(x, sval)

dt = .5 * delta_eta1

for i = 1:nc_eta1
   xi = eta1_min + (i-1)*delta_eta1
   ex_exact[i] = - cos_k(xi)*dt
end

compute_rhs_from_function!(rho, maxwell_1d, cos_k, deg-1)
fill!(ex, 0.0)
compute_e_from_j!(ex, maxwell_1d, dt .* rho, 1 )
sval =  eval_uniform_periodic_spline_curve(deg-1, ex);

plot(x, ex_exact)
scatter!(x, sval)

# +
mode     = 2
eta1_min = .0
eta1_max = 2π
nc_eta1  = 256

Lx = eta1_max - eta1_min

delta_eta1 = Lx / nc_eta1

domain = [eta1_min, eta1_max]

deg = 3

maxwell_1d = Maxwell1DFEM(domain, nc_eta1, deg)

ex = zeros(Float64, nc_eta1)
ey = zeros(Float64, nc_eta1)
bz = zeros(Float64, nc_eta1)

bz_exact = zeros(Float64, nc_eta1)
ex_exact = zeros(Float64, nc_eta1)
ey_exact = zeros(Float64, nc_eta1)
rho      = zeros(Float64, nc_eta1)
sval     = zeros(Float64, nc_eta1)

cos_k(x) = cos(mode*2*pi*x/Lx) 
sin_k(x) = sin(mode*2*pi*x/Lx) 


# Test Maxwell on By and Ez 
#--------------------------
# Set time stepping parameters

time  = 0.0
dt    = .5 * delta_eta1
nstep = 10

# Compute initial fields 
ex = 0.0 # 0-form -> splines of degree deg
l2projection!(bz, maxwell_1d, cos_k, deg-1) # 0-form -> splines of degree deg-1

plot(bz)

# +
#for istep = 1:nstep 

   compute_b_from_e!( bz, maxwell_1d, 0.5*dt, ey)

plot(bz)
# -

compute_e_from_b!( ey, maxwell_1d,     dt, bz)
plot(ey)

# +
compute_b_from_e!( bz, maxwell_1d, 0.5*dt, ey)

time = time + dt

for i = 1:nc_eta1
   xi = eta1_min + (i-1)*delta_eta1
   ey_exact[i] = sin(mode * 2π * xi/Lx) * sin(mode * 2π * time/Lx)
   bz_exact[i] = cos(mode * 2π * xi/Lx) * cos(mode * 2π * time/Lx)
end

sval = eval_uniform_periodic_spline_curve(deg, ey)
@show err_ey = norm(sval .- ey_exact)

sval = eval_uniform_periodic_spline_curve(deg-1, bz)
@show err_bz = norm(sval .- bz_exact)

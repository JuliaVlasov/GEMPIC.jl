
@testset " Maxwell 1D FEM solver " begin

using GEMPIC
using LinearAlgebra

"""
Test 1D Maxwell spline finite element solver on a periodic grid

`L_x` domain dimensions and M is an integer.

```math
B_z(x,y,t) = \\cos(\\frac{2 M \\pi}{L_x} x)  \\cos(\\frac{2 M \\pi}{L_x} t)
```
```math
E_y(x,y,t) = \\sin(\\frac{2 M \\pi}{L_x} x)  \\sin(\\frac{2 M \\pi}{L_x} t)
```

"""

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

# Test Poisson
#-------------
# Set exact solution
for i = 1:nc_eta1
   xi = eta1_min + (i-1)*delta_eta1
   ex_exact[i] = sin_k(xi)/(2.0*mode*pi/Lx)
end

compute_rhs_from_function!( rho, maxwell_1d, cos_k, deg)

compute_e_from_rho!( ex, maxwell_1d, rho ) 

# Evaluate spline curve at grid points and compute error
# Ex is a 1-form, i.e. one spline degree lower
sval = eval_uniform_periodic_spline_curve(deg-1, ex)
err_ex = maximum(abs.(sval .- ex_exact))
println( " error Poisson  $err_ex ")
@test err_ex ≈ 0.0 atol = 1e-6

# Test Ampere
#-------------
# Set time step
dt = .5 * delta_eta1
# Set exact solution
for i = 1:nc_eta1
   xi = eta1_min + (i-1)*delta_eta1
   ex_exact[i] = - cos_k(xi)*dt
end

compute_rhs_from_function!(rho, maxwell_1d, cos_k, deg-1)
fill!(ex, 0.0)
compute_e_from_j!(ex, maxwell_1d, dt .* rho, 1 )

# Evaluate spline curve at grid points and compute error
# Ex is a 1-form, i.e. one spline degree lower
sval =  eval_uniform_periodic_spline_curve(deg-1, ex)
err_ex2 = maximum(abs.(sval .- ex_exact))
#println( " error Ampere  $err_ex2 ")
@test err_ex2 ≈ 0.0 atol = 1e-6

l2norm =  l2norm_squared(maxwell_1d, ex, deg-1)
err_l2norm = l2norm - dt^2 * pi
#println( " error l2 norm $err_l2norm ")

# Test Maxwell on By and Ez 
#--------------------------
# Set time stepping parameters

time  = 0.0
dt    = .5 * delta_eta1
nstep = 10

# Compute initial fields 
ex = 0.0 # 0-form -> splines of degree deg
l2projection!(bz, maxwell_1d, cos_k, deg-1) # 0-form -> splines of degree deg-1

for istep = 1:nstep 

   compute_b_from_e!( bz, maxwell_1d, 0.5*dt, ey)
   compute_e_from_b!( ey, maxwell_1d,     dt, bz)
   compute_b_from_e!( bz, maxwell_1d, 0.5*dt, ey)
   
   time = time + dt

   for i = 1:nc_eta1
      xi = eta1_min + (i-1)*delta_eta1
      ey_exact[i] = sin(mode * 2π * xi/Lx) * sin(mode * 2π * time/Lx)
      bz_exact[i] = cos(mode * 2π * xi/Lx) * cos(mode * 2π * time/Lx)
   end

   sval = eval_uniform_periodic_spline_curve(deg, ey)
   err_ey = norm(sval .- ey_exact)

   sval = eval_uniform_periodic_spline_curve(deg-1, bz)
   err_bz = norm(sval .- bz_exact)

   @test err_ey ≈ 0.0 atol = 1e-2
   @test err_bz ≈ 0.0 atol = 1e-2

end

end

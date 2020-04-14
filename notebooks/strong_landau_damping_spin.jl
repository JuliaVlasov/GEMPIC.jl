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
#     display_name: Julia 1.4.0
#     language: julia
#     name: julia-1.4
# ---

# # Strong Landau Damping
#
# Electrostatic example of strong Landau damping
# $$
# f(x,v) =\frac{1}{2\pi\sigma^2}  \exp 
# \Big( - \frac{v_1^2 + v_2^2}{2\sigma^2} \Big)
# ( 1+\alpha \cos(k x)􏰁),
# $$
#

using ProgressMeter, Plots
using CSV, Dates, FFTW
using GEMPIC

# +
function run( steps :: Int64) 
   
   σ, μ = 0.02, 0.0
   kx, α = 1.004355, 0.001
   xmin, xmax = 0, 2π/kx
   nx = 1024 
   n_particles = 100000
   mesh = Mesh( xmin, xmax, nx)
   spline_degree = 3
   
   df = CosSumGaussianSpin([[kx]], [α], [[σ]], [[μ]] )
   
   particle_group2 = ParticleGroup{1,1}( n_particles, n_spin=3)   
   sampler = ParticleSampler{1,1}( :sobol, n_particles)
   
   sample!(particle_group2, sampler, df, mesh)
   
   particle_group = ParticleGroup{1,1}( n_particles, n_spin=3)   
   GEMPIC.set_common_weight(particle_group, (1.0/n_particles))

   for  i_part = 1:n_particles

       x = zeros( 1 )
       v = zeros( 1 )
       s = zeros( 3 )
       w = zeros( 1 )
       x = GEMPIC.get_x(particle_group2, i_part)
       v = GEMPIC.get_v(particle_group2, i_part)

       s1 = GEMPIC.get_spin(particle_group2, i_part, 1)
       s2 = GEMPIC.get_spin(particle_group2, i_part, 2)
       s3 = GEMPIC.get_spin(particle_group2, i_part, 3)

       w = GEMPIC.get_weights(particle_group2, i_part)

       GEMPIC.set_x(particle_group, i_part, x[1])
       GEMPIC.set_v(particle_group, i_part, v[1])
       GEMPIC.set_spin(particle_group, i_part, 1, s1)
       GEMPIC.set_spin(particle_group, i_part, 2, s2)
       GEMPIC.set_spin(particle_group, i_part, 3, s3)
       GEMPIC.set_weights(particle_group, i_part, w[1])

   end
   
   kernel_smoother2 = ParticleMeshCoupling( mesh, n_particles, spline_degree-2, :galerkin) 
   kernel_smoother1 = ParticleMeshCoupling( mesh, n_particles, spline_degree-1, :galerkin)    
   kernel_smoother0 = ParticleMeshCoupling( mesh, n_particles, spline_degree, :galerkin)
   
   rho = zeros(Float64, nx)
   efield_poisson = zeros(Float64, nx)
   
   # Init!ialize the field solver
   maxwell_solver = Maxwell1DFEM(mesh, spline_degree)
   # efield by Poisson
   solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
   
   # Initialize the arrays for the spline coefficients of the fields
   k0 = 12.0523 
   E0 = 10
   ww = 12.104827940833333

   Ey(x) = E0*cos(k0*x)
   Ez(x) = E0*sin(k0*x)
   Ay(x) = -E0/ww*sin(k0*x)
   Az(x) = E0/ww*cos(k0*x)
   
   efield_dofs = [efield_poisson, zeros(Float64, nx), zeros(Float64, nx)]
   afield_dofs = [zeros(Float64, nx), zeros(Float64, nx)]
   
   l2projection!( efield_dofs[2], maxwell_solver, Ey, spline_degree)
   l2projection!( efield_dofs[3], maxwell_solver, Ez, spline_degree)
   l2projection!( afield_dofs[1], maxwell_solver, Ay, spline_degree)
   l2projection!( afield_dofs[2], maxwell_solver, Az, spline_degree)
       
   propagator = HamiltonianSplittingSpin( maxwell_solver,
                                      kernel_smoother0, 
                                      kernel_smoother1, 
                                      kernel_smoother2,
                                      particle_group,
                                      efield_dofs,
                                      afield_dofs)
  
   efield_dofs_n = propagator.e_dofs
   
   thdiag = TimeHistoryDiagnosticsSpin( particle_group, maxwell_solver, 
                           kernel_smoother0, kernel_smoother1 );
   
   Δt = 0.002

   @showprogress 1 for jstep = 1:steps # loop over time
   
       # Strang splitting
       strang_splitting!(propagator, Δt, 1)
   
       solve_poisson!( efield_poisson, particle_group, 
                       kernel_smoother0, maxwell_solver, rho)
       
       write_step!(thdiag, jstep * Δt, spline_degree, 
                       efield_dofs,  afield_dofs,
                       efield_dofs_n, efield_poisson, propagator)

       if (jstep % 1000 == 0)
           GEMPIC.save( "particles", jstep, particle_group)
       end
       
   end

   thdiag.data

end
# +
results = run(10000) # choose number of steps

CSV.write("thdiag-$(now()).csv", results)
# -

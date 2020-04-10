# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.1
#   kernelspec:
#     display_name: Julia 1.3.0
#     language: julia
#     name: julia-1.3
# ---

# +
using ProgressMeter, Plots
using CSV, Dates
using GEMPIC

function run( steps) 
   
   # # Strong Landau Damping
   #
   # Electrostatic example of strong Landau damping
   # $$
   # f(x,v) =\frac{1}{2\pi\sigma^2}  \exp 
   # \Big( - \frac{v_1^2 + v_2^2}{2\sigma^2} \Big)
   # ( 1+\alpha \cos(k x)􏰁),
   # $$
   #
   # The physical parameters 
   
   # +
   σ, μ = 0.02, 0.0
   kx, α = 1.004355, 0.001
   xmin, xmax = 0, 2π/kx
   domain = [xmin, xmax, xmax - xmin]
   nx = 1024 
   n_particles = 100000
   mesh = Mesh( xmin, xmax, nx)
   spline_degree = 3
   
   df = SpinCosSumGaussian{1,1,3}([[kx]], [α], [[σ]], [[μ]] )
   
   mass, charge = 1.0, 1.0
   
   particle_group2 = SpinParticleGroup{1,1,3}( n_particles, mass, charge, 1)   
   sampler = SpinParticleSampler{1,1,3}( :sobol, n_particles)
   
   sample!(particle_group2, sampler, df, mesh)
   
   particle_group = SpinParticleGroup{1,1,3}( n_particles, mass, charge, 1)   
   set_common_weight(particle_group, (1.0/n_particles))
   for  i_part = 1:n_particles
       x = zeros( 1 )
       v = zeros( 1 )
       s = zeros( 3 )
       w = zeros( 1 )
       x = get_x(particle_group2, i_part)
       v = get_v(particle_group2, i_part)
       s[1] =  get_s1(particle_group2, i_part)
       s[2] =  get_s2(particle_group2, i_part)
       s[3] =  get_s3(particle_group2, i_part)
       w = get_weights(particle_group2, i_part)
       set_x(particle_group, i_part, x[1])
       set_v(particle_group, i_part, v[1])
       set_s1(particle_group, i_part, s[1])
       set_s2(particle_group, i_part, s[2])
       set_s3(particle_group, i_part, s[3])
       set_weights(particle_group, i_part, w[1])
   end
   
   kernel_smoother2 = SpinParticleMeshCoupling( domain, [nx], n_particles, spline_degree-2, :galerkin) 
   kernel_smoother1 = SpinParticleMeshCoupling( domain, [nx], n_particles, spline_degree-1, :galerkin)    
   kernel_smoother0 = SpinParticleMeshCoupling( domain, [nx], n_particles, spline_degree, :galerkin)
   
   efield_poisson = zeros(Float64, nx)
   
   # Init!ialize the field solver
   maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
   # efield by Poisson
   solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
   
   # +
   # # +
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
       
   propagator = SpinHamiltonianSplitting( maxwell_solver,
                                      kernel_smoother0, 
                                      kernel_smoother1, 
                                      kernel_smoother2,
                                      particle_group,
                                      efield_dofs,
                                      afield_dofs,
                                      domain, nx);
   
   efield_dofs_n = propagator.e_dofs
   
   thdiag = SpinTimeHistoryDiagnostics( particle_group, maxwell_solver, 
                           kernel_smoother0, kernel_smoother1 );
   
   Δt = 0.002

   mode1 = zeros(ComplexF64,steps)
   mode2 = zeros(ComplexF64,steps)
   store = zeros(ComplexF64,nx)
   
   ss11 = Float64[]
   ss12 = Float64[]
   ss13 = Float64[]
   ss21 = Float64[]
   ss22 = Float64[]
   ss23 = Float64[]
   ss31 = Float64[]
   ss32 = Float64[]
   ss33 = Float64[]

   electric = zeros(ComplexF64, steps, nx)
   
   @showprogress 1 for j = 1:steps # loop over time
   
       # Strang splitting
       strang_splitting!(propagator, Δt, 1)
   
       solve_poisson!( efield_poisson, particle_group, 
                       kernel_smoother0, maxwell_solver, rho)
       
       write_step!(thdiag, j * Δt, spline_degree, 
                       efield_dofs,  afield_dofs,
                       efield_dofs_n, efield_poisson, propagator)

       for i in eachindex(store)
           xi = 2π/kx/nx*(i-1)
           store[i] = evaluate(propagator.kernel_smoother_1, xi, propagator.e_dofs[1])
       end

       fft!(store)
       electric[j,:] .= store 
       
       push!(ss11, get_s1(propagator.particle_group, 1)) 
       push!(ss12, get_s2(propagator.particle_group, 1))
       push!(ss13, get_s3(propagator.particle_group, 1)) 

       push!(ss21, get_s1(propagator.particle_group, 100)) 
       push!(ss22, get_s2(propagator.particle_group, 100)) 
       push!(ss23, get_s3(propagator.particle_group, 100)) 

       push!(ss31, get_s1(propagator.particle_group, 800)) 
       push!(ss32, get_s2(propagator.particle_group, 800)) 
       push!(ss33, get_s3(propagator.particle_group, 800)) 
       
   end

   return thdiag.data

end

results = run(10000)

CSV.write("thdiag-$(now()).csv", results)

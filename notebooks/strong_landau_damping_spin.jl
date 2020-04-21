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
using Sobol
using Random
using Distributions
using JLD2
using FileIO


function main_sample!( pg   :: ParticleGroup{1,1}, 
                  ps   :: ParticleSampler, 
                  df   :: CosSumGaussianSpin, 
                  mesh :: Mesh )

    s = zeros( pg.n_spin )

    theta = 0.0
    phi = 0.0
    n_rnds = 0
    if df.params.n_gaussians > 1
       n_rnds = 1
    end

    δ = zeros(df.params.n_gaussians)
    for i_v=1:df.params.n_gaussians
       δ[i_v] = sum(df.params.δ[1:i_v])
    end
    
    n_rnds += 2
    rdn = zeros(3)
    
    if ps.sampling_type == :sobol
       rng_sobol  = SobolSeq(1)
    end 

    rng_random = MersenneTwister(ps.seed)

    d = Normal()

    for i_part = 1:(pg.n_particles)  
       
       if ps.sampling_type == :sobol
           x = mesh.xmin[1] + Sobol.next!(rng_sobol)[1] * mesh.Lx[1]
       else
           x = mesh.xmin[1] + rand(rng_random) * mesh.Lx[1]
       end

       v = rand(rng_random, d)

       # For multiple Gaussian, draw which one to take
       rnd_no = rdn[3]
        
       i_gauss = 1
       while( rnd_no > δ[i_gauss] )
          i_gauss += 1
       end
       v = v * df.params.σ[i_gauss] + df.params.μ[i_gauss]

#for a peaked initial condition in the s direction 
       s = [0, 0, 1]
#for a uniformly distributed initial condition on the sphere
#       for tt = 1:10
#            s[1] = randn()
#            s[2] = randn()
#            s[3] = randn()
#            if norm(s) > 10^(-4)
#                break
#            end
#        end
#        s .= s./norm(s)

       # Set weight according to value of perturbation
       w  = GEMPIC.eval_x_density(df, x) * prod(mesh.Lx) 
        
       # Copy the generated numbers to the particle
       GEMPIC.set_x(pg, i_part, x[1])
       GEMPIC.set_v(pg, i_part, v)
       GEMPIC.set_spin(pg, i_part, 1, s[1])
       GEMPIC.set_spin(pg, i_part, 2, s[2])
       GEMPIC.set_spin(pg, i_part, 3, s[3])
       # Set weights.
       GEMPIC.set_weights(pg, i_part, w)
        
    end
       
end



# +
function run( steps :: Int64) 
   
   σ, μ = 0.02, 0.0
   kx, α = 1.004355, 0.001
   xmin, xmax = 0, 2π/kx
   nx = 32 #NC 1024 
   n_particles = 1000 #NC 100000
   mesh = Mesh( xmin, xmax, nx)
   spline_degree = 3
   
   df = CosSumGaussianSpin([[kx]], [α], [[σ]], [[μ]] )
   
   particle_group2 = ParticleGroup{1,1}( n_particles, n_spin=3)   
   sampler = ParticleSampler{1,1}( :sobol, n_particles)
   
   main_sample!(particle_group2, sampler, df, mesh)
   
   particle_group = ParticleGroup{1,1}( n_particles, n_spin=3)
   
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

   #Test 1 with only Ey and Ay (and with HH=0)
   Ey(x) = E0*cos(k0*x)
   Ez(x) = 0.0*E0*sin(k0*x)
   Ay(x) = -E0/ww*sin(k0*x)
   Az(x) = 0.0*E0/ww*cos(k0*x)
   
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

#    mode = zeros(ComplexF64,steps,nx)
     th_modes = Vector{ComplexF64}[]
    elec_tmp = zeros(Float64,nx)
   
   Δt = 0.002

   @showprogress 1 for jstep = 1:steps # loop over time
   
       # Strang splitting
       strang_splitting!(propagator, Δt, 1)
   
       solve_poisson!( efield_poisson, particle_group, 
                       kernel_smoother0, maxwell_solver, rho)
       
       write_step!(thdiag, jstep * Δt, spline_degree, 
                       efield_dofs,  afield_dofs,
                       efield_dofs_n, efield_poisson, propagator)

       #diagnostics

       #store particles at some specific times 
       if (jstep % 1000 == 0)
           GEMPIC.save( "save_particles", jstep, particle_group)
       end

       #Fourier modes of the longitudinal electric field
       for i = 1:nx
           elec_tmp[i] = GEMPIC.evaluate(thdiag.kernel_smoother_1, (i-1)*propagator.delta_x[1], efield_dofs[1])
       end
       push!(th_modes,fft(elec_tmp))


   end

   thdiag.data, th_modes

end
# +
thdiag, th_modes = run(100) # choose number of steps

CSV.write("thdiag-$(now()).csv", thdiag)

file="th_modes"
#datafile = @sprintf("%s.jld2", file)

FileIO.save("th_modes.jld2", Dict("modes" => th_modes))


# -




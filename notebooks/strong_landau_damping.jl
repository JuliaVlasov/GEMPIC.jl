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

using ProgressMeter, Plots
using GEMPIC

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
σ, μ = 1.0, 0.0
kx, α = 0.5, 0.5
xmin, xmax = 0, 2π/kx
∆t = 0.05
nx = 32 
n_particles = 100000
mesh = Mesh( xmin, xmax, nx)
spline_degree = 3

df = CosSumGaussian{1,2}([[kx]], [α], [[σ,σ]], [[μ,μ]] )

mass, charge = 1.0, 1.0
particle_group = ParticleGroup{1,2}(n_particles)   
sampler = ParticleSampler{1,2}( :sobol, true, n_particles)

sample!(particle_group, sampler, df, mesh)

xp = view(particle_group.array, 1, :)
vp = view(particle_group.array, 2:3, :)
wp = view(particle_group.array, 4, :);
# -

p = plot(layout=(3,1))
histogram!(p[1,1], xp, weights=wp, normalize=true, bins = 100, lab = "")
plot!(p[1,1], x-> (1+α*cos(kx*x))/(2π/kx), 0., 2π/kx, lab="")
histogram!(p[2,1], vp[1,:], weights=wp, normalize=true, bins = 100, lab = "")
plot!(p[2,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
histogram!(p[3,1], vp[2,:], weights=wp, normalize=true, bins = 100, lab = "")
plot!(p[3,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
savefig("histograms.png")


kernel_smoother1 = ParticleMeshCoupling( mesh, n_particles, spline_degree-1, :galerkin)    
kernel_smoother0 = ParticleMeshCoupling( mesh, n_particles, spline_degree, :galerkin)
rho = zeros(Float64, nx)
xg = LinRange(xmin, xmax, nx)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, rho)
plot( xg, sval )
savefig("charge_density.png")

efield_poisson = zeros(Float64, nx)
# Initialize the field solver
maxwell_solver = Maxwell1DFEM(mesh, spline_degree)
# efield by Poisson
solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
sval = eval_uniform_periodic_spline_curve(spline_degree-1, efield_poisson)
plot( xg, sval )       
savefig("electric_field.png")

# Initialize the arrays for the spline coefficients of the fields

# +
function run( steps)

    σ, μ = 1.0, 0.0
    kx, α = 0.5, 0.5
    xmin, xmax = 0, 2π/kx
    ∆t = 0.05
    nx = 32 
    n_particles = 100000
    mesh = Mesh( xmin, xmax, nx)
    spline_degree = 3
    
    df = CosSumGaussian{1,2}([[kx]], [α], [[σ,σ]], [[μ,μ]] )
    
    mass, charge = 1.0, 1.0
    particle_group = ParticleGroup{1,2}(n_particles)   
    sampler = ParticleSampler{1,2}( :sobol, true, n_particles)
    
    sample!(particle_group, sampler, df, mesh)
    
    kernel_smoother1 = ParticleMeshCoupling( mesh, n_particles, spline_degree-1, :galerkin)    
    kernel_smoother0 = ParticleMeshCoupling( mesh, n_particles, spline_degree, :galerkin)
    
    rho = zeros(Float64, nx)
    efield_poisson = zeros(Float64, nx)
    # Initialize the field solver
    maxwell_solver = Maxwell1DFEM(mesh, spline_degree)
    # efield by Poisson
    solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
    
    efield_dofs = [efield_poisson, zeros(Float64, nx)]
    bfield_dofs = zeros(Float64, nx)
        
    propagator = HamiltonianSplitting( maxwell_solver,
                                       kernel_smoother0, 
                                       kernel_smoother1, 
                                       particle_group,
                                       efield_dofs,
                                       bfield_dofs)
    
    efield_dofs_n = propagator.e_dofs
    
    thdiag = TimeHistoryDiagnostics( particle_group, maxwell_solver, 
                            kernel_smoother0, kernel_smoother1 )
    
    Δt = 0.01
    
    @showprogress 1 for j = 1:steps # loop over time
    
        # Strang splitting
        strang_splitting!(propagator, Δt, 1)
    
        # Diagnostics
        solve_poisson!( efield_poisson, particle_group, 
                        kernel_smoother0, maxwell_solver, rho)
        
        write_step!(thdiag, j * Δt, spline_degree, 
                        efield_dofs,  bfield_dofs,
                        efield_dofs_n, efield_poisson)
    
    end
    
    return thdiag.data
    
end
# -

@time results = run(5000) # change number of steps

plot(results[!,:Time], log.(results[!,:PotentialEnergyE1]))

savefig("potential_energy.png")



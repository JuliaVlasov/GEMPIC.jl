# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Julia 1.1.1
#     language: julia
#     name: julia-1.1
# ---

# +
using ProgressMeter, Plots

include("../src/mesh.jl")
include("../src/distributions.jl")
include("../src/low_level_bsplines.jl")
include("../src/splinepp.jl")
include("../src/maxwell_1d_fem.jl")
include("../src/particle_group.jl")
include("../src/particle_mesh_coupling.jl")
include("../src/hamiltonian_splitting.jl")
include("../src/particle_sampling.jl")
include("../src/diagnostics.jl")

# -

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
k, α = 0.5, 1.0
xmin, xmax = 0, 2π/k
domain = [xmin, xmax, xmax - xmin]
∆t = 0.05
nx = 32 
n_particles = 100000
mesh = Mesh( xmin, xmax, nx)
spline_degree = 3

df = CosSumGaussian{1,2}([[k]], [α], [[σ,σ]], [[μ,μ]] )


mass, charge = 1.0, 1.0
particle_group = ParticleGroup{1,2}( n_particles, mass, charge, 1)   
sampler = ParticleSampler{1,2}( :sobol, true, n_particles)

sample!(particle_group, sampler, df, mesh)

xp = Vector{Float64}[] # particles data
for i in 1:n_particles
    push!(xp, vcat(get_x(particle_group,i), 
            get_v(particle_group,i),
            get_weights(particle_group,i)))
end

# +
p = vcat(xp[1:end]'...);

scatter(p[:,2], p[:,3], markersize=1)

p = histogram([ p[:,1],p[:,2]], normalize=true, bins = 100,  layout=(2,1), lab = "draws")
    
# -

kernel_smoother1 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree-1, :galerkin)    
kernel_smoother0 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree, :galerkin)
rho = zeros(Float64, nx)
ex = zeros(Float64, nx)
maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
solve_poisson!( ex, particle_group, 
                kernel_smoother0, maxwell_solver, rho)
xg = LinRange(xmin, xmax, nx)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, rho)
plot( xg, sval )

efield_poisson = zeros(Float64, nx)
# Init!ialize the field solver
maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
# efield by Poisson
solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
sval = eval_uniform_periodic_spline_curve(spline_degree-1, efield_poisson)
plot( xg, sval )       

# +
# # +
# Initialize the arrays for the spline coefficients of the fields
efield_dofs = [efield_poisson, zeros(Float64, nx)]
bfield_dofs = zeros(Float64, nx)
    
propagator = HamiltonianSplitting( maxwell_solver,
                                   kernel_smoother0, 
                                   kernel_smoother1, 
                                   particle_group,
                                   efield_dofs,
                                   bfield_dofs,
                                   domain);

efield_dofs_n = propagator.e_dofs

thdiag = TimeHistoryDiagnostics( particle_group, maxwell_solver, 
                        kernel_smoother0, kernel_smoother1 );

steps, Δt = 500, 0.05

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

# +
first(thdiag.data, 10)

using Gadfly


Gadfly.plot(thdiag.data, x=:Time, y=:PotentialEnergyE2, Geom.point, Geom.line)
# -



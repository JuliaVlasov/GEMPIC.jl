using BenchmarkTools

import Base.Threads: @sync, @spawn, nthreads, threadid

include("../src/mesh.jl")
include("../src/low_level_bsplines.jl")
include("../src/splinepp.jl")
include("../src/distributions.jl")
include("../src/maxwell_1d_fem.jl")
include("../src/particle_group.jl")
include("../src/particle_mesh_coupling.jl")
include("../src/particle_sampling.jl")
include("../src/landau_damping.jl")
include("../src/hamiltonian_splitting.jl")
include("../src/hamiltonian_splitting_boris.jl")
include("../src/diagnostics.jl")

const β = 0.0001
const k = 1.25
const α = 0.0
const σ = [0.2,  0.005773502691896]
const μ = [0.0, 0.0]

const nx   = 64
const xmin = 0.0
const xmax = 2π / k
const n_particles    = 100000
const sampling_case  = :sobol
const symmetric      = true
const splitting_case = :symplectic
const spline_degree  = 3

function setup( )

    mesh   = Mesh( xmin, xmax, nx)
    
    beta_sin_k(x) = beta * sin(2π * x / (xmax - xmin)) 
    
    df = CosSumGaussian{1,2}([[k]], [α], [σ], [μ] )
    
     # Initialize the particles   (mass and charge set to 1.0 with one weight)
    particle_group = ParticleGroup{1,2}( n_particles)   
    
    sampler = ParticleSampler{1,2}( sampling_case, symmetric, n_particles)
    
    sample!(  particle_group, sampler, df, mesh)
    
    kernel_smoother1 = ParticleMeshCoupling( mesh, n_particles, spline_degree-1, :galerkin)    
    kernel_smoother0 = ParticleMeshCoupling( mesh, n_particles, spline_degree, :galerkin)
    rho = zeros(Float64, nx);
    efield_poisson = zeros(Float64, nx)
    # Initialize the field solver
    maxwell_solver = Maxwell1DFEM(mesh, spline_degree)
    # efield by Poisson
    solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )    
    
    # Initialize the arrays for the spline coefficients of the fields
    efield_dofs = [copy(efield_poisson), zeros(Float64, nx)]
    bfield_dofs = zeros(Float64, nx)
        
    propagator = HamiltonianSplitting( maxwell_solver,
                                       kernel_smoother0, 
                                       kernel_smoother1, 
                                       particle_group,
                                       efield_dofs,
                                       bfield_dofs)
    
    efield_dofs_n = propagator.e_dofs
    
    beta_cos_k(x) = β * cos(2π * x / (xmax - xmin))
    l2projection!( bfield_dofs, maxwell_solver, beta_cos_k, spline_degree-1)
    
    thdiag = TimeHistoryDiagnostics( particle_group, maxwell_solver, 
                            kernel_smoother0, kernel_smoother1 );
    
    return propagator
    
end

propagator = setup()

strang_splitting!(propagator, 0.05, 1) # trigger compilation
@time  strang_splitting!(propagator, 0.05, 10)

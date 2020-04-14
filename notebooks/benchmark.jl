using GEMPIC

import Base.Threads: @sync, @spawn, nthreads, threadid

function setup( )

    β = 0.0001
    k = 1.25
    α = 0.0
    σ = [0.2,  0.005773502691896]
    μ = [0.0, 0.0]
    
    nx   = 32
    xmin = 0.0
    xmax = 2π / k
    n_particles    = 100000
    sampling_case  = :sobol
    symmetric      = true
    splitting_case = :symplectic
    spline_degree  = 3

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

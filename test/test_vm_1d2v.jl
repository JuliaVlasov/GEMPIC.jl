import GEMPIC: get_charge, get_x, add_charge!, get_mass
import GEMPIC: inner_product, evaluate

"""
Simulation of 1d2v Vlasov-Maxwell with simple PIC method, 
periodic boundary conditions, Weibel instability. 
FEM with splines, degree 3 for B and 2 for E
"""

@testset " PIC VM 1D2V " begin

    delta_t         = 0.05
    n_time_steps    = 5
    beta            = 0.0001
    initial_bfield  = :cos
    
    k = [[1.25]]
    α = [0.0]
    σ = [[0.2, 0.005773502691896]]
    μ = [[0.0, 0.0]]
    
    nx   = 32
    xmin = 0.0
    xmax = 5.02654824574
    
    n_particles    = 100000
    sampling_case  = :sobol
    symmetric      = true
    splitting_case = :symplectic
    spline_degree  = 3
    
    # TODO mesh and domain contain same information
    mesh   = Mesh( xmin, xmax, nx)
    domain = [xmin, xmax, xmax - xmin ]

    beta_cos_k(x) = beta * cos(2π * x / domain[3]) 
    beta_sin_k(x) = beta * sin(2π * x / domain[3]) 
    
    degree_smoother = spline_degree

    # Distribution function 
    df = CosSumGaussian{1,2}( k, α, σ, μ )
    
    sampler = ParticleSampler{1,2}( sampling_case, symmetric, n_particles)
    
    # Initialize the particles : mass and charge set to 1.0 with one weight
    particle_group = ParticleGroup{1,2}( n_particles, 1.0, 1.0, 1)
    
    # Initialize the field solver
    maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
    
    kernel_smoother1 = ParticleMeshCoupling( domain, [nx], n_particles, 
                 spline_degree, :galerkin)
    
    kernel_smoother0 = ParticleMeshCoupling( domain, [nx], n_particles, 
                 spline_degree, :galerkin)
    
    # Initialize the arrays for the spline coefficients of the fields
    efield_dofs = [zeros(Float64, nx), zeros(Float64, nx)]
    bfield_dofs = zeros(Float64, nx)

    # Set the initial fields
    rho = zeros(Float64, nx)
    efield_poisson = zeros(Float64, nx)
    
    propagator = HamiltonianSplitting( maxwell_solver,
                                       kernel_smoother0, 
                                       kernel_smoother1, 
                                       particle_group,
                                       efield_dofs, 
                                       bfield_dofs,
                                       domain )

    efield_dofs_n = propagator.e_dofs

    sample!( particle_group, sampler, df, mesh )

    # efield 1 by Poisson, rho is computed inside the function
    solve_poisson!( efield_dofs[1], particle_group, kernel_smoother0, 
                    maxwell_solver, rho )

    # bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))

    if initial_bfield == :cos
        l2projection!( bfield_dofs, maxwell_solver, beta_cos_k, degree_smoother-1)
    else
        l2projection!( bfield_dofs, maxwell_solver, beta_sin_k, degree_smoother-1)
    end
    
    thdiag = TimeHistoryDiagnostics( particle_group, maxwell_solver,
                                     kernel_smoother0, kernel_smoother1 )

    write_step!( thdiag, 0.0, degree_smoother, 
                 efield_dofs, bfield_dofs, efield_dofs_n, efield_poisson)

    for j = 1:n_time_steps # loop over time

       # Strang splitting
       strang_splitting!(propagator, delta_t, 1)

       # Diagnostics
       solve_poisson!( efield_poisson, particle_group, 
                       kernel_smoother0, maxwell_solver, rho)

       write_step!( thdiag, j * delta_t, degree_smoother, 
                    efield_dofs, bfield_dofs, efield_dofs, efield_poisson)

       @test true

    end

    # Compute final rho
    fill!(rho, 0.0)
    for i_part = 1:particle_group.n_particles
       xi = get_x(particle_group, i_part)
       wi = get_charge( particle_group, i_part)
       add_charge!(rho, kernel_smoother0, xi, wi)
    end

    @show thdiag.data

end

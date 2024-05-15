using ProgressMeter, Plots
using CSV, DataFrames
using GEMPIC

function run( steps)

    σ, μ = 1.0, 0.0
    kx, α = 0.5, 0.5
    xmin, xmax = 0, 2π/kx
    dt = 0.05
    nx = 32 
    n_particles = 100000
    mesh = OneDGrid( xmin, xmax, nx)
    spline_degree = 3
    
    mass, charge = 1.0, 1.0
    particle_group = ParticleGroup{1,1}(n_particles)   

    sample!(particle_group, α, kx, σ, mesh)
    
    kernel_smoother1 = ParticleMeshCoupling1D( mesh, n_particles, spline_degree-1, :galerkin)    
    kernel_smoother0 = ParticleMeshCoupling1D( mesh, n_particles, spline_degree, :galerkin)
    
    rho = zeros(Float64, nx)
    efield_poisson = zeros(Float64, nx)
    # Initialize the field solver
    maxwell_solver = Maxwell1DFEM(mesh, spline_degree)
    # efield by Poisson
    solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
    
    efield_dofs = [efield_poisson, zeros(Float64, nx)]
    bfield_dofs = zeros(Float64, nx)
        
    propagator = HamiltonianSplitting{1,1}( maxwell_solver,
                                       kernel_smoother0, 
                                       kernel_smoother1, 
                                       particle_group,
                                       efield_dofs,
                                       bfield_dofs)
    
    efield_dofs_n = propagator.e_dofs
    
    thdiag = TimeHistoryDiagnostics( particle_group, maxwell_solver, 
                            kernel_smoother0, kernel_smoother1 )
    
    @showprogress 1 for j = 1:steps # loop over time
    
        # Strang splitting
        strang_splitting!(propagator, dt, 1)
    
        # Diagnostics
        solve_poisson!( efield_poisson, particle_group, 
                        kernel_smoother0, maxwell_solver, rho)
        
        write_step!(thdiag, j * dt, spline_degree, 
                        efield_dofs,  bfield_dofs,
                        efield_dofs_n, efield_poisson)
    
    end
    
    return thdiag.data
    
end

@time results = run(1000) # change number of steps
CSV.write("results.csv", results)

plot(results.Time, log.(results.PotentialEnergyE1))

png("results")

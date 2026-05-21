using GEMPIC
using Plots
using ProgressMeter

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

    time = Float64[]
    energy = Float64[]
    
    @showprogress 1 for j = 1:steps # loop over time
    
        # Strang splitting
        strang_splitting!(propagator, dt, 1)
    
        # Diagnostics
        solve_poisson!( efield_poisson, particle_group, 
                        kernel_smoother0, maxwell_solver, rho)

        push!(time, j * dt)
        push!(energy, GEMPIC.inner_product(maxwell_solver, efield_dofs[1], efield_dofs_n[1], spline_degree - 1))

        
    end
    
    return time, energy
    
end

@time time, energy = run(1000) # change number of steps

plot(time, log.(energy))

png("results")

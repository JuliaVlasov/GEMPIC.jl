@testset "VM 1D1V" begin

    σ, μ = 1.0, 0.0
    kx, α = 0.5, 0.5
    xmin, xmax = 0, 2π / kx
    n_particles = 100_000
    nx = 32
    mesh = OneDGrid(xmin, xmax, nx)
    spline_degree = 3

    df = CosSumGaussian{1,1}([[kx]], [α], [[σ]], [[μ]])

    particle_group = ParticleGroup{1,1}(n_particles)
    sampler = ParticleSampler{1,1}(:sobol, n_particles)

    sample!(particle_group, sampler, df, mesh)

    kernel_smoother0 = ParticleMeshCoupling1D(mesh, n_particles, spline_degree, :galerkin)
    kernel_smoother1 = ParticleMeshCoupling1D(mesh, n_particles, spline_degree - 1, :galerkin)

    rho = zeros(Float64, nx)
    efield_poisson = zeros(Float64, nx)

    # Initialize the field solver
    maxwell_solver = Maxwell1DFEM(mesh, spline_degree)
    # efield by Poisson
    solve_poisson!(efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho)

    efield_dofs = [zeros(Float64, nx) for _ in 1:3]
    bfield_dofs = zeros(Float64, nx)

    efield_dofs[1] .= efield_poisson

    propagator = HamiltonianSplitting{1,1}(
        maxwell_solver,
        kernel_smoother0,
        kernel_smoother1,
        particle_group,
        efield_dofs, 
        bfield_dofs
    )

    efield_dofs_n = propagator.e_dofs

    thdiag = TimeHistoryDiagnostics(
        particle_group, maxwell_solver, kernel_smoother0, kernel_smoother1
    )

    steps, Δt = 2, 0.1

    # Strang splitting
    strang_splitting!(propagator, Δt, 1)

    solve_poisson!(efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho)

    strang_splitting!(propagator, Δt, 1)

    solve_poisson!(efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho)

    @test true

end

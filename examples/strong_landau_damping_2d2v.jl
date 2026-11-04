using GEMPIC
using Plots

function run_strong_landau_damping_2d2v(;
    steps=500, n_particles=180_000, nx=32, ny=32, dt=0.02
)
    xmin, xmax = 0.0, 4π
    ymin, ymax = 0.0, 4π
    degree = 3
    mesh = TwoDGrid(xmin, xmax, nx, ymin, ymax, ny)
    distribution = CosSumGaussian{2,2}(
        [[0.5, 0.0], [0.0, 0.5]],
        [0.5, 0.5],
        [[1.0, 1.0]],
        [[0.0, 0.0]],
    )

    particles = ParticleGroup{2,2}(n_particles)
    sampler = ParticleSampler{2,2}(:sobol, n_particles)
    sample!(particles, sampler, distribution, mesh)

    maxwell = TwoDMaxwell(mesh, degree)
    density_kernel = ParticleMeshCoupling2D(
        particles, mesh, degree, :galerkin
    )
    rho = zeros(nx * ny)
    efield = [zeros(nx * ny) for _ in 1:3]
    bfield = [zeros(nx * ny) for _ in 1:3]
    solve_poisson!(efield, particles, density_kernel, maxwell, rho)

    splitting = HamiltonianSplitting2D(maxwell, particles, efield, bfield)
    diagnostics = TimeHistoryDiagnostics2D(particles, maxwell)
    for step in 1:steps
        strang_splitting!(splitting, dt, 1)
        write_step!(diagnostics, step * dt, efield, bfield)
    end
    return diagnostics.data
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = run_strong_landau_damping_2d2v()
    show(results, allrows=true)
    println()
    electric_energy = results.ElectricEnergyEx .+ results.ElectricEnergyEy
    energy_plot = plot(
        results.Time,
        log.(electric_energy);
        xlabel="Time",
        ylabel="log Electric Energy",
        label="Ex + Ey",
        title="2d2v Strong Landau Damping",
    )
    savefig(energy_plot, "strong_landau_damping_2d2v_electric_energy.png")
end

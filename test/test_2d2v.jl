@testset "2d2v simulation" begin
    mesh = TwoDGrid(0.0, 2π, 8, 0.0, 2π, 8)
    distribution = CosSumGaussian{2,2}(
        [[1.0, 0.0], [0.0, 1.0]],
        [0.05, 0.05],
        [[1.0, 1.0]],
        [[0.0, 0.0]],
    )
    sampled = ParticleGroup{2,2}(4096)
    sampler = ParticleSampler{2,2}(:sobol, 4096)
    sample!(sampled, sampler, distribution, mesh)
    @test all((sampled.array[1, :] .>= mesh.xmin) .& (sampled.array[1, :] .< mesh.xmax))
    @test all((sampled.array[2, :] .>= mesh.ymin) .& (sampled.array[2, :] .< mesh.ymax))
    @test all(isfinite, sampled.array)
    degree = 2
    maxwell = TwoDMaxwell(mesh, degree)
    density_kernel = ParticleMeshCoupling2D(sampled, mesh, degree, :galerkin)
    rho = zeros(mesh.nx * mesh.ny)
    efield = [zeros(mesh.nx * mesh.ny) for _ in 1:3]
    solve_poisson!(efield, sampled, density_kernel, maxwell, rho)
    @test sum(rho) ≈ mesh.dimx * mesh.dimy rtol = 2e-3
    initial_electric_energy =
        GEMPIC.inner_product(maxwell, efield[1], efield[1], 1, 1) +
        GEMPIC.inner_product(maxwell, efield[2], efield[2], 2, 1)
    @test initial_electric_energy ≈ 0.05^2 * mesh.dimx * mesh.dimy rtol = 0.05

    particles = ParticleGroup{2,2}(5; common_weight=0.2)
    for i in 1:particles.n_particles
        particles.array[:, i] .= [2π - 0.01, π, 0.25, -0.1, 1.0]
    end
    density = ParticleMeshCoupling2D(particles, mesh, degree, :galerkin)
    rho = zeros(mesh.nx * mesh.ny)
    charge_kernel = GEMPIC.ParticleMeshCoupling2DAnisotropic(
        mesh, degree, degree, :galerkin
    )
    GEMPIC.add_charge!(rho, charge_kernel, 0.2, 0.3, 1.0)
    @test sum(rho) ≈ 1.0

    efield = [zeros(mesh.nx * mesh.ny) for _ in 1:3]
    bfield = [zeros(mesh.nx * mesh.ny) for _ in 1:3]
    solve_poisson!(efield, particles, density, maxwell, rho)
    splitting = HamiltonianSplitting2D(maxwell, particles, efield, bfield)

    x0 = particles.array[1, 1]
    GEMPIC.operatorHp1(splitting, 0.1)
    @test particles.array[1, 1] ≈ mod(x0 + 0.1 * 0.25, 2π)
    @test sum(splitting.current_dofs[1]) ≈ 0.025

    strang_splitting!(splitting, 0.01, 2)
    @test all(isfinite, particles.array)
    @test all(isfinite, efield[1])
    @test all(isfinite, efield[2])
    @test all(iszero, efield[3])
    @test all(iszero, bfield[1])
    @test all(iszero, bfield[2])

    diagnostics = TimeHistoryDiagnostics2D(particles, maxwell)
    write_step!(diagnostics, 0.02, efield, bfield)
    @test size(diagnostics.data, 1) == 1
    @test all(isfinite, collect(diagnostics.data[1, :]))
    @test diagnostics.data.GaussLawError[1] < 1e-10
end

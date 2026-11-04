@testset "1d1v diagnostics" begin
    import GEMPIC: get_v, get_mass
    nx, np, degree = 16, 64, 3
    mesh = OneDGrid(0.0, 4π, nx)
    pg = ParticleGroup{1,1}(np; charge=-2.0, mass=3.0, common_weight=0.25)
    for i in 1:np
        pg.array[1, i] = 4π * (i - 0.5) / np
        pg.array[2, i] = sin(i / 7)
        pg.array[3, i] = 1.0 + 0.2cos(i / 5)
    end
    kernel0 = ParticleMeshCoupling1D(mesh, np, degree, :galerkin)
    kernel1 = ParticleMeshCoupling1D(mesh, np, degree - 1, :galerkin)
    solver = Maxwell1DFEM(mesh, degree)
    rho, e1 = zeros(nx), zeros(nx)
    solve_poisson!(e1, pg, kernel0, solver, rho)
    efields = [copy(e1), zeros(nx)]
    bfield = zeros(nx)
    diagnostics = TimeHistoryDiagnostics(pg, solver, kernel0, kernel1)

    transfer = sum(
        pg.array[2, i] * GEMPIC.evaluate(kernel1, pg.array[1, i], efields[1]) *
        (pg.array[3, i] * pg.charge * pg.common_weight) for i in 1:np
    )
    @test GEMPIC.pic_diagnostics_transfer(pg, kernel0, kernel1, efields) ≈ transfer
    @test GEMPIC.pic_diagnostics_vvb(pg, kernel1, bfield) == 0.0

    @test write_step!(diagnostics, 0.0, degree, efields, bfield, efields, e1) !== nothing
    row = diagnostics.data[1, :]
    expected_kinetic = sum(
        pg.array[2, i]^2 * pg.array[3, i] * pg.mass * pg.common_weight for i in 1:np
    )
    expected_momentum = sum(
        pg.array[2, i] * pg.array[3, i] * pg.mass * pg.common_weight for i in 1:np
    )
    @test row.KineticEnergy ≈ expected_kinetic
    @test row.Momentum1 ≈ expected_momentum
    @test row.Momentum2 == 0.0
    @test row.PotentialEnergyE1 > 0.0
    @test row.PotentialEnergyE2 == 0.0
    @test row.PotentialEnergyB3 == 0.0
    @test row.Transfer ≈ transfer
    @test row.VVB == 0.0
    @test row.Poynting == 0.0
    @test row.ErrorPoisson == 0.0
    expected_component_kinetic = sum(
        get_v(pg, i)[1]^2 * get_mass(pg, i)[1] for i in 1:np
    )
    @test GEMPIC.pic_diagnostics_hpi(pg, 1, 0.0) ≈ expected_component_kinetic
    @test names(diagnostics.data) == [
        "Time", "KineticEnergy", "Momentum1", "Momentum2", "PotentialEnergyE1",
        "PotentialEnergyE2", "PotentialEnergyB3", "Transfer", "VVB", "Poynting",
        "ErrorPoisson",
    ]
end

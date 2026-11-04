function reference_charge!(rho, pg, kernel)
    fill!(rho, 0.0)
    for i in 1:pg.n_particles
        GEMPIC.add_charge!(rho, kernel, pg.array[1, i], GEMPIC.get_charge(pg, i))
    end
    return rho
end

function charge_deposition_bytes(rho, pg, kernel)
    GEMPIC.deposit_charge!(rho, pg, kernel)
    return @allocated GEMPIC.deposit_charge!(rho, pg, kernel)
end

@testset "Cell-moment charge deposition" begin
    nx = 16
    mesh = OneDGrid(-2.0, 6.0, nx)
    for V in (1, 2), degree in 1:5, smoothing in (:galerkin, :collocation)
        np = 1024
        pg = ParticleGroup{1,V}(np; charge=-1.5, mass=2.0, common_weight=0.25, n_weights=2)
        for i in 1:np
            pg.array[1, i] = -2.0 + 8.0 * (i - 1) / (np - 1)
            pg.array[V + 2, i] = 1.0 + 0.5sin(i)
            pg.array[V + 3, i] = 100.0
        end
        kernel = ParticleMeshCoupling1D(mesh, np, degree, smoothing)
        rho, reference = zeros(nx), zeros(nx)
        moments = kernel.charge_moments
        for layout in (:uniform, :clustered, :signed, :periodic)
            if layout == :clustered
                pg.array[1, :] .= -1.73
            elseif layout == :signed
                pg.array[V + 2, :] .= sin.(1:np)
            elseif layout == :periodic
                pg.array[1, :] .= range(-10.0, 14.0; length=np)
            end
            reference_charge!(reference, pg, kernel)
            fill!(moments, NaN)
            @test GEMPIC.deposit_charge!(rho, pg, kernel) === rho
            @test rho ≈ reference rtol=2e-13 atol=2e-12
            expected_charge = pg.charge * pg.common_weight * sum(pg.array[V + 2, :])
            @test sum(rho) ≈ kernel.scaling * expected_charge rtol=2e-13 atol=2e-12
            @test kernel.charge_moments === moments
        end
        @test charge_deposition_bytes(rho, pg, kernel) == 0
    end

    pg = ParticleGroup{1,2}(0)
    kernel = ParticleMeshCoupling1D(mesh, 0, 3, :galerkin)
    rho = ones(nx)
    GEMPIC.deposit_charge!(rho, pg, kernel)
    @test all(iszero, rho)
end

@testset "Poisson bulk deposition" begin
    mesh = OneDGrid(-2.0, 6.0, 32)
    for V in (1, 2), degree in (2, 3)
        np = 2048
        pg = ParticleGroup{1,V}(np; charge=-1.5, common_weight=0.25)
        for i in 1:np
            pg.array[1, i] = -2.0 + 8.0 * (i - 1) / (np - 1)
            pg.array[V + 2, i] = 1.0 + 0.5cos(i)
        end
        kernel = ParticleMeshCoupling1D(mesh, np, degree, :galerkin)
        solver = Maxwell1DFEM(mesh, degree)
        rho, reference = zeros(32), zeros(32)
        e, e_reference = zeros(32), zeros(32)
        original_particles = copy(pg.array)
        for _ in 1:3
            reference_charge!(reference, pg, kernel)
            GEMPIC.compute_e_from_rho!(e_reference, solver, reference)
            solve_poisson!(e, pg, kernel, solver, rho)
            @test rho ≈ reference rtol=2e-13 atol=2e-12
            @test e ≈ e_reference rtol=2e-13 atol=2e-12
            @test pg.array == original_particles
        end
        @test (@allocated solve_poisson!(e, pg, kernel, solver, rho)) == 0
    end
end

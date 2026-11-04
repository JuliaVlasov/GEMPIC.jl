function serial_he!(h, dt)
    pg = h.particle_group
    for i in 1:pg.n_particles
        x = pg.array[1, i]
        e1 = GEMPIC.evaluate(h.kernel_smoother_1, x, h.e_dofs[1])
        e2 = GEMPIC.evaluate(h.kernel_smoother_0, x, h.e_dofs[2])
        pg.array[2, i] += dt * pg.q_over_m * e1
        pg.array[3, i] += dt * pg.q_over_m * e2
    end
    GEMPIC.compute_b_from_e!(h.b_dofs, h.maxwell_solver, dt, h.e_dofs[2])
end

function he_conversion_bytes(pp, spline, b)
    GEMPIC.b_to_pp!(pp, spline, length(b), b)
    return @allocated GEMPIC.b_to_pp!(pp, spline, length(b), b)
end

@testset "Periodic pp field evaluation" begin
    nx = 16
    mesh = OneDGrid(-2.0, 6.0, nx)
    b = sin.(1:nx) .+ cos.(2 .* (1:nx))
    for degree in 1:5
        kernel = ParticleMeshCoupling1D(mesh, 1, degree, :galerkin)
        pp = fill(NaN, degree + 1, nx)
        @test GEMPIC.b_to_pp!(pp, kernel.spline_pp, nx, b) === pp
        for cell in 0:(nx - 1), offset in (0.0, 0.125, 0.5, 0.875)
            x = mesh.xmin + (cell + offset) * kernel.delta_x
            @test GEMPIC.evaluate_pp(kernel, x, pp) ≈ GEMPIC.evaluate(kernel, x, b) atol=1e-14
        end
        for x in (-10.125, -2.0, prevfloat(6.0), 6.0, 6.125, 14.0)
            @test GEMPIC.evaluate_pp(kernel, x, pp) ≈ GEMPIC.evaluate(kernel, x, b) atol=1e-14
        end
        @test he_conversion_bytes(pp, kernel.spline_pp, b) == 0
    end
end

@testset "HE read-only polynomial workspace" begin
    nx = 32
    np = 256 * Threads.nthreads()
    mesh = OneDGrid(-2.0, 6.0, nx)
    for degree in (2, 3)
        pg = ParticleGroup{1,2}(np; charge=-1.0, mass=2.0)
        for i in 1:np
            pg.array[:, i] .= (-2.0 + 8.0 * (i - 1) / (np - 1), sin(i), cos(i), 1.0)
        end
        h = HamiltonianSplitting{1,2}(
            Maxwell1DFEM(mesh, degree),
            ParticleMeshCoupling1D(mesh, np, degree, :galerkin),
            ParticleMeshCoupling1D(mesh, np, degree - 1, :galerkin),
            pg, [sin.(1:nx), cos.(1:nx)], sin.(2 .* (1:nx)),
        )
        reference = deepcopy(h)
        buffers = h.he_pp
        original_positions = copy(pg.array[1, :])
        original_weights = copy(pg.array[4, :])
        for dt in (0.0, 0.025, -0.025, 0.4)
            # Refresh coefficients after field changes, not just at construction.
            h.e_dofs[1] .+= 0.1
            h.e_dofs[2] .-= 0.2
            reference.e_dofs[1] .= h.e_dofs[1]
            reference.e_dofs[2] .= h.e_dofs[2]
            fill!.(buffers, NaN)
            for kernel in (h.kernel_smoother_0, h.kernel_smoother_1)
                fill!(kernel.spline_val, NaN)
                fill!(kernel.spline_val_more, NaN)
            end
            serial_he!(reference, dt)
            GEMPIC.operatorHE(h, dt)
            @test pg.array ≈ reference.particle_group.array atol=1e-13
            @test pg.array[1, :] == original_positions
            @test pg.array[4, :] == original_weights
            @test h.b_dofs ≈ reference.b_dofs atol=1e-14
            @test h.e_dofs == reference.e_dofs
            @test h.he_pp[1] === buffers[1]
            @test h.he_pp[2] === buffers[2]
            @test all(isnan, h.kernel_smoother_0.spline_val)
            @test all(isnan, h.kernel_smoother_1.spline_val)
        end
        chunk = first(h.chunks)
        GEMPIC.push_he!(h, 0.025, chunk)
        @test (@allocated GEMPIC.push_he!(h, 0.025, chunk)) == 0
        if Threads.nthreads() == 1
            GEMPIC.operatorHE(h, 0.025)
            @test (@allocated GEMPIC.operatorHE(h, 0.025)) == 0
        end
    end
end

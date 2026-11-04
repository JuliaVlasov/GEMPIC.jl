function serial_hp1!(h, dt)
    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)
    pg = h.particle_group
    for i in 1:pg.n_particles
        x = pg.array[1, i]
        x_new = x + dt * pg.array[2, i]
        w = pg.array[4, i] * pg.charge * pg.common_weight
        pg.array[3, i] = GEMPIC.add_current_update_v!(
            h.j_dofs[1], h.kernel_smoother_1, x, x_new, w, pg.q_over_m,
            h.b_dofs, pg.array[3, i],
        )
        pg.array[1, i] = mod(x_new, h.Lx)
    end
    for i in 1:pg.n_particles
        w = pg.array[4, i] * pg.charge * pg.common_weight
        GEMPIC.add_charge!(h.j_dofs[2], h.kernel_smoother_0, pg.array[1, i], w)
    end
    GEMPIC.compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.j_dofs[1], 1)
end

@testset "Hp1 reusable chunk workspace" begin
    np = 256 * Threads.nthreads()
    nx = 32
    mesh = OneDGrid(0.0, 4π, nx)
    for degree in (2, 3)
        pg = ParticleGroup{1,2}(np; charge=-1.0, mass=2.0)
        for i in 1:np
            pg.array[:, i] .= (4π * (i - 0.5) / np, 3sin(i), cos(i), 1 + 0.5cos(i))
        end
        h = HamiltonianSplitting{1,2}(
            Maxwell1DFEM(mesh, degree),
            ParticleMeshCoupling1D(mesh, np, degree, :galerkin),
            ParticleMeshCoupling1D(mesh, np, degree - 1, :galerkin),
            pg, [sin.(1:nx), cos.(1:nx)], sin.(1:nx),
        )
        reference = deepcopy(h)
        workspaces = copy(h.hp1_workspaces)
        for workspace in workspaces
            @test workspace.kernel_0.spline_val !== h.kernel_smoother_0.spline_val
            @test workspace.kernel_1.spline_val !== h.kernel_smoother_1.spline_val
        end
        if length(workspaces) > 1
            @test workspaces[1].kernel_1.spline_val !== workspaces[2].kernel_1.spline_val
        end

        # Exercise zero motion, both directions, cell crossings and periodic wrapping.
        for dt in (0.0, 0.05, 0.7, -0.7, 0.05)
            for workspace in workspaces
                fill!(workspace.current, NaN)
                fill!(workspace.charge, NaN)
            end
            serial_hp1!(reference, dt)
            GEMPIC.operatorHp1(h, dt)
            @test h.particle_group.array ≈ reference.particle_group.array atol=1e-13
            @test h.j_dofs[1] ≈ reference.j_dofs[1] atol=1e-13
            @test h.j_dofs[2] ≈ reference.j_dofs[2] atol=1e-13
            @test h.e_dofs[1] ≈ reference.e_dofs[1] atol=1e-13
            @test h.e_dofs[2] == reference.e_dofs[2]
            @test all(a === b for (a, b) in zip(workspaces, h.hp1_workspaces))
        end

        workspace = first(workspaces)
        chunk = first(h.chunks)
        GEMPIC.deposit_hp1!(h, 0.05, chunk, workspace)
        @test (@allocated GEMPIC.deposit_hp1!(h, 0.05, chunk, workspace)) == 0
    end
end

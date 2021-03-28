@testset "ParticleMeshCoupling 2D" begin
    n_cells = 10 # Number of cells
    n_particles = 4 # Number of particles
    spline_degree = 3 # Spline degree
    xmin, xmax, nx = 0.0, 2.0, n_cells
    ymin, ymax, ny = 0.0, 1.0, n_cells

    mesh = TwoDGrid(xmin, xmax, nx, ymin, ymax, ny)

    volume = (xmax - xmin) * (ymax - ymin)

    x_vec = [0.1 0.65 0.7 1.5; 0.0 0.0 0.0 0.0]
    v_vec = [1.5 0.00 0.0 0.0; 0.0 0.5 0.0 0.0]

    pg = ParticleGroup{2,2}(n_particles; charge=1.0, mass=1.0, n_weights=1)

    for i_part in 1:n_particles
        pg.array[1, i_part] = x_vec[1, i_part]
        pg.array[2, i_part] = x_vec[2, i_part]
        pg.array[3, i_part] = v_vec[1, i_part]
        pg.array[4, i_part] = v_vec[2, i_part]
        pg.array[5, i_part] = 1.0 / n_particles
    end

    # Initialize the kernel
    kernel = ParticleMeshCoupling2D(pg, mesh, spline_degree, :collocation)

    # Reference values of the shape factors
    index_grid = zeros(Int, (2, 4))
    index_grid[1, :] .= [-2, 1, 1, 5]
    index_grid[2, :] .= [-3, -3, -3, -3]

    v_grid = zeros(Float64, (4, 2, 4))
    v_grid[:, 1, 1] .= [
        2.0833333333333332E-002,
        0.47916666666666663,
        0.47916666666666663,
        2.0833333333333332E-002,
    ]
    v_grid[:, 1, 3] .= v_grid[:, 1, 1]
    v_grid[:, 1, 4] .= v_grid[:, 1, 1]
    v_grid[:, 1, 2] .= [
        7.0312500000000000E-002,
        0.61197916666666663,
        0.31510416666666663,
        2.6041666666666665E-003,
    ]
    v_grid[1, 2, :] .= 0.0
    v_grid[2, 2, :] .= 1.0 / 6.0
    v_grid[3, 2, :] .= 2.0 / 3.0
    v_grid[4, 2, :] .= 1.0 / 6.0

    ρ_dofs = zeros(Float64, 100)
    ρ_dofs1 = zeros(Float64, 100)
    ρ_dofs_ref = zeros(Float64, 100)
    ρ_dofs_pp = zeros(Float64, (16, 100))

    # Accumulate ρ
    for i_part in 1:n_particles
        x = pg.array[1, i_part]
        y = pg.array[2, i_part]
        w = pg.array[5, i_part]
        add_charge!(ρ_dofs, kernel, x, y, w)
        add_charge_pp!(ρ_dofs1, kernel, x, y, w)
    end

    ρ_dofs_ref[8:10] .= v_grid[1:3, 1, 1]
    ρ_dofs_ref[1] = v_grid[4, 1, 1]
    ρ_dofs_ref[1:4] .= ρ_dofs_ref[1:4] .+ v_grid[:, 1, 2] .+ v_grid[:, 1, 3]
    ρ_dofs_ref[5:8] .= ρ_dofs_ref[5:8] .+ v_grid[:, 1, 4]

    ρ_dofs_ref[71:80] .= ρ_dofs_ref[1:10] ./ 6.0
    ρ_dofs_ref[81:90] .= ρ_dofs_ref[1:10] .* 2.0 / 3.0
    ρ_dofs_ref[91:100] .= ρ_dofs_ref[1:10] ./ 6.0
    ρ_dofs_ref[1:10] .= 0.0

    ρ_dofs_ref .*= n_cells^2 / volume / n_particles

    @test ρ_dofs ≈ ρ_dofs_ref
    @test ρ_dofs1 ≈ ρ_dofs_ref

    GEMPIC.b_to_pp_2d!(ρ_dofs_pp, kernel.spline1, kernel.spline2, ρ_dofs)

    # For evaluation check
    particle_values = zeros(4)
    particle_values1 = zeros(4)
    particle_values_ref = zeros(4)

    # Test function evaluation
    for i_part in 1:n_particles
        xp = pg.array[1, i_part]
        yp = pg.array[2, i_part]
        particle_values[i_part] = evaluate(kernel, xp, yp, ρ_dofs)
        particle_values1[i_part] = evaluate_pp(kernel, xp, yp, ρ_dofs_pp)
    end

    particle_values_ref =
        [1.1560058593749998, 2.3149278428819446, 2.2656250000000000, 1.1512586805555554] ./
        volume

    fill!(particle_values_ref, 0.0)

    for i_part in 1:n_particles
        for i in 1:4
            i1 = mod(index_grid[1, i_part] + i - 2, n_cells)
            for j in 1:4
                i2 = mod(index_grid[2, i_part] + j - 2, n_cells)
                res =
                    v_grid[i, 1, i_part] *
                    v_grid[j, 2, i_part] *
                    ρ_dofs_ref[i1 + i2 * n_cells + 1]
                particle_values_ref[i_part] += res
            end
        end
    end

    @test particle_values ≈ particle_values_ref
    @test particle_values1 ≈ particle_values_ref
end

using Test
using GEMPIC

@testset "Maxwell 2D" begin
    function evaluate_spline_2d(nx1, nx2, degs, dofs)
        deg1, deg2 = degs
        vals = collect(reshape(dofs, nx1, nx2))
        for j in 1:nx2
            vals[:, j] .= GEMPIC.eval_uniform_periodic_spline_curve(deg1, vals[:, j])
        end
        for i in 1:nx1
            vals[i, :] .= GEMPIC.eval_uniform_periodic_spline_curve(deg2, vals[i, :])
        end
        return vec(vals)
    end

    x1min, x1max = 0.0, 2π
    nx1 = 16
    x2min, x2max = 0.0, 2π
    nx2 = 32

    mesh = TwoDGrid(x1min, x1max, nx1, x2min, x2max, nx2)

    deg = 3
    delta_t = 0.01
    nsteps = 300

    maxwell = TwoDMaxwell(mesh, deg)

    efield = [zeros(nx1 * nx2) for _ in 1:3]
    bfield = deepcopy(efield)
    efield_val = deepcopy(efield)
    bfield_val = deepcopy(efield)
    efield_ref = deepcopy(efield)
    bfield_ref = deepcopy(efield)

    x = LinRange(x1min, x1max, nx1 + 1)[1:(end - 1)] .* transpose(ones(nx2))
    y = ones(nx1) .* transpose(LinRange(x2min, x2max, nx2 + 1)[1:(end - 1)])

    w1 = sqrt(3)
    w2 = sqrt(3)

    rho = zeros(nx1 * nx2)
    rho_ref = zeros(nx1 * nx2)
    time = 0.0

    sin_k = (x, y) -> sin((x + y) - w1 * time)
    cos_k = (x, y) -> cos((x + y) - w1 * time)

    rho .= compute_rhs_from_function(maxwell, cos_k, 1, 0)

    compute_e_from_rho!(efield, maxwell, rho)

    efield_val1 = evaluate_spline_2d(nx1, nx2, (deg - 1, deg), efield[1])
    efield_val2 = evaluate_spline_2d(nx1, nx2, (deg, deg - 1), efield[2])
    efield_val3 = evaluate_spline_2d(nx1, nx2, (deg, deg), efield[3])

    efield_ref[1] .= vec(sin_k.(x, y) ./ 2)
    efield_ref[2] .= efield_ref[1]
    efield_ref[3] .= 0.0

    @test efield_val1 ≈ efield_ref[1] rtol = 1e-4
    @test efield_val2 ≈ efield_ref[2] rtol = 1e-4
    @test efield_val3 ≈ efield_ref[3]

    e1(x, y) = cos(x) * sin(y) * sin(sqrt(2) * time) / sqrt(2)
    e2(x, y) = -sin(x) * cos(y) * sin(sqrt(2) * time) / sqrt(2)
    b3(x, y) = -cos(x) * cos(y) * cos(sqrt(2) * time)

    time = -0.5 * delta_t

    bfield[1] .= l2projection(maxwell, e1, 1, 2)
    bfield[2] .= l2projection(maxwell, e2, 2, 2)
    bfield[3] .= l2projection(maxwell, b3, 3, 2)

    time = 0.0

    efield[1] .= l2projection(maxwell, e1, 1, 1)
    efield[2] .= l2projection(maxwell, e2, 2, 1)
    efield[3] .= l2projection(maxwell, b3, 3, 1)
    efield[3] .*= -1

    for istep in 1:nsteps
        compute_b_from_e!(bfield, maxwell, delta_t, efield)
        compute_e_from_b!(efield, maxwell, delta_t, bfield)
    end

    bfield_val[1] = evaluate_spline_2d(nx1, nx2, (deg, deg - 1), bfield[1])
    bfield_val[2] = evaluate_spline_2d(nx1, nx2, (deg - 1, deg), bfield[2])
    bfield_val[3] = evaluate_spline_2d(nx1, nx2, (deg - 1, deg - 1), bfield[3])

    time = (nsteps - 0.5) * delta_t
    bfield_ref[1] .= vec(e1.(x, y))
    bfield_ref[2] .= vec(e2.(x, y))
    bfield_ref[3] .= vec(b3.(x, y))

    @test bfield_ref[1] ≈ bfield_val[1] rtol = 1e-4
    @test bfield_ref[2] ≈ bfield_val[2] rtol = 1e-4
    @test bfield_ref[3] ≈ bfield_val[3] rtol = 1e-3

    efield_val[1] = evaluate_spline_2d(nx1, nx2, (deg - 1, deg), efield[1])
    efield_val[2] = evaluate_spline_2d(nx1, nx2, (deg, deg - 1), efield[2])
    efield_val[3] = evaluate_spline_2d(nx1, nx2, (deg, deg), efield[3])

    time = nsteps * delta_t

    efield_ref[1] .= vec(e1.(x, y))
    efield_ref[2] .= vec(e2.(x, y))
    efield_ref[3] .= -vec(b3.(x, y))

    @test efield_ref[1] ≈ efield_val[1] rtol = 1e-4
    @test efield_ref[2] ≈ efield_val[2] rtol = 1e-4
    @test efield_ref[3] ≈ efield_val[3] rtol = 1e-3

    efield[1] = l2projection(maxwell, cos_k, 1, 1)

    error2 = GEMPIC.inner_product(maxwell, efield[1], efield[1], 1, 1) - 2 * pi^2
    println(" Error in L2 norm squared: $error2")
    @test error2 ≈ 0 atol = 1e-5

    rho .= compute_rhs_from_function(maxwell, sin_k, 1, 1)

    compute_e_from_j!(efield[1], maxwell, rho, 1)

    efield_val1 = evaluate_spline_2d(nx1, nx2, [deg - 1, deg], efield[1])

    for i in eachindex(x, y)
        efield_ref[1][i] = cos_k(x[i], y[i]) - sin_k(x[i], y[i])
    end

    @show maximum(abs.(efield_val1 .- efield_ref[1]))
    @test efield_val1 ≈ efield_ref[1] atol = 1e-2

    time = 0.0
    rho_ref = compute_rhs_from_function(maxwell, cos_k, 1, 0)
    rho_ref .*= 2.0

    efield[1] = l2projection(maxwell, sin_k, 1, 1)
    efield[2] = l2projection(maxwell, sin_k, 2, 1)
    efield[3] = l2projection(maxwell, sin_k, 3, 1)

    compute_rho_from_e!(rho, maxwell, efield)

    @show maximum(abs.(rho .- rho_ref))

    @test rho ≈ rho_ref
end

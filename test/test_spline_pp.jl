@testset " spline pp 1d " begin
    import GEMPIC: SplinePP, b_to_pp, uniform_bsplines_eval_basis, horner_1d

    ncells = 8

    domain = [0.0, 2π]
    delta_x = (domain[2] - domain[1]) / ncells

    b = rand(ncells)

    for d in 1:5
        xp = 4.8141437173987800

        spline_pp = SplinePP(d, ncells)

        pp = b_to_pp(spline_pp, ncells, b)

        xi = (xp - domain[1]) / delta_x
        index = floor(Int64, xi) + 1
        xi = xi - (index - 1)
        res1 = horner_1d(d, pp, xi, index)

        index = index - d

        val = uniform_bsplines_eval_basis(d, xi)

        res2 = 0.0

        for i in 1:(d + 1)
            index1d = (index + i - 2) % ncells + 1
            res2 += b[index1d] * val[i]
        end

        @test res1 ≈ res2

        @test spline_pp.degree ≈ d

        xp = rand()
        index = 1
        res1 = horner_1d(d, pp, xp, index)

        res2 = 0.0
        for i in 1:(d + 1)
            res2 += pp[i, 1] * xp^((d + 1) - i)
        end

        @test res1 ≈ res2
    end
end

@testset " spline pp 2d " begin
    using Random

    ncells = 50
    d = 3

    b = zeros(ncells * ncells)
    pp = zeros((d + 1) * (d + 1), ncells * ncells)

    dx = 2π / ncells
    dy = 2π / ncells

    rng = MersenneTwister(42)
    rand!(rng, b)

    spline1 = SplinePP(d, ncells)
    spline2 = SplinePP(d, ncells)

    GEMPIC.b_to_pp_2d!(pp, spline1, spline2, b)

    xp, yp = rand(rng, 2) .* 2π

    xi = xp / dx
    yi = yp / dy

    ind_x = floor(Int, xi) + 1
    ind_y = floor(Int, yi) + 1

    xi -= (ind_x - 1)
    yi -= (ind_y - 1)

    res1 = GEMPIC.horner_2d((d, d), pp, (xi, yi), (ind_x, ind_y), (ncells, ncells))

    ind_x -= d
    ind_y -= d

    val1 = uniform_bsplines_eval_basis(d, xi)
    val2 = uniform_bsplines_eval_basis(d, yi)

    res2 = 0.0
    for i in 1:(d + 1)
        idx1 = mod(ind_x + i - 2, ncells)
        for j in 1:(d + 1)
            idx2 = mod(ind_y + j - 2, ncells)
            index2d = idx1 + idx2 * ncells + 1
            res2 += b[index2d] * val1[i] * val2[j]
        end
    end

    @test res1 ≈ res2

    xp = rand(rng, 2)

    res1 = GEMPIC.horner_2d((d, d), pp, xp, [1, 1], [1, 1])

    res2 = 0.0
    for i in 1:(d + 1), j in 1:(d + 1)
        res2 += pp[i + (j - 1) * (d + 1), 1] * xp[1]^((d + 1) - i) * xp[2]^((d + 1) - j)
    end

    @test res1 ≈ res2
end

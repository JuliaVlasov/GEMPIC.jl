@testset "Spin Sampling" begin

    function test_sampling_spin(
        sampling_type::Symbol,
        pg::ParticleGroup{D,V},
        df::CosSumGaussianSpin,
    ) where {D,V}

        mean = zeros(2)
        sigma = zeros(2)

        n_particles = pg.n_particles

        sampling = ParticleSampler{D,V}(sampling_type, n_particles)

        sample!(pg, sampling, df, mesh)

        for i_part = 1:n_particles
            xi = get_x(pg, i_part)
            vi = get_v(pg, i_part)
            mean[1] += xi[1]
            mean[2] += vi[1]
        end

        mean = mean / n_particles

        for i_part = 1:n_particles
            xi = get_x(pg, i_part)
            vi = get_v(pg, i_part)
            sigma[1] += (xi[1] - mean[1])^2
            sigma[2] += (vi[1] - mean[2])^2
        end

        sigma = sigma / (n_particles - 1)

        mean, sigma

    end

    σ, μ = 0.02, 0.0
    kx, α = 1.004355, 0.001
    n_particles = 100000
    xmin = 0.0
    xmax = 2π / kx
    Lx = xmax - xmin
    nx = 512

    mesh = OneDGrid(xmax, xmin, nx)

    pg = ParticleGroup{1,1}(n_particles, n_spin = 3)

    params = (k = [[kx]], α = [α], σ = [[σ]], μ = [[μ]])

    df = CosSumGaussianSpin(params...)

    mean, sigma = test_sampling_spin(:sobol, pg, df)

    @test true

end

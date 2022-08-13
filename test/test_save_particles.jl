@testset "Save particles data into file" begin
    n_particles = 100000
    xmin = 1.0::Float64
    xmax = 4π + 1.0
    Lx = xmax - xmin
    nx = 64

    mesh = OneDGrid(xmax, xmin, nx)

    pg = ParticleGroup{1,2}(n_particles)

    params = (k=[[0.5]], α=[0.01], σ=[[0.1, 2.0]], μ=[[0.0, 0.0]])

    df = CosSumGaussian{1,2}(params...)

    n_particles = pg.n_particles

    sampler = ParticleSampler{1,2}(:sobol, true, n_particles)

    sample!(pg, sampler, df, mesh)

    GEMPIC.save("test", 1, pg)

    @test isfile("test-000001.jld2")

end

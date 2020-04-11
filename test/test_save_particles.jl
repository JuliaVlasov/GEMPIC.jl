@testset "Save particles data into file" begin

n_particles = 100000
xmin        = 1.0 :: Float64
xmax        = 4π + 1.0
Lx          = xmax - xmin  
nx          = 64

mesh = Mesh( xmax, xmin, nx)

pg = ParticleGroup{1,2}(n_particles, 1.0, 1.0, 1)

params = ( k = [[0.5]],
           α = [0.01],
           σ = [[0.1, 2.0]],
           μ = [[0.0, 0.0]]
)

df = CosSumGaussian{1,2}(params...)

n_particles = pg.n_particles
   
sampler = ParticleSampler{1, 2}( :sobol, true, n_particles )

sample!( pg, sampler, df, mesh )

GEMPIC.save( "test", 1, pg)

@test isfile("test-000001.jld2")

σ, μ = 0.02, 0.0
kx, α = 1.004355, 0.001
xmin, xmax = 0, 2π/kx
domain = [xmin, xmax, xmax - xmin]
nx = 64 
n_particles = 1000
mesh = Mesh( xmin, xmax, nx)
spline_degree = 3
   
df = SpinCosSumGaussian{1,1,3}([[kx]], [α], [[σ]], [[μ]] )
   
mass, charge = 1.0, 1.0
   
spg = SpinParticleGroup{1,1,3}( n_particles, mass, charge, 1)   
sampler = SpinParticleSampler{1,1,3}( :sobol, n_particles)
   
sample!(spg, sampler, df, mesh)
   
GEMPIC.save( "test", 2, spg)

@test isfile("test-000002.jld2")

end

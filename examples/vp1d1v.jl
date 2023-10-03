using ParticleInCell
using Plots
using Random

dt = 0.1
nsteps = 100
alpha = 0.5
kx = 0.5

nx = 64
xmin, xmax = 0.0, 2π / kx

n_particles = 1000000
degree_smoother = 3

mesh = OneDGrid( xmin, xmax, nx)

particles = ParticleGroup{1,1}( n_particles, charge=1.0, mass=1.0, n_weights=1)

sampler = LandauDamping( alpha, kx )

ParticleInCell.sample!( particles, mesh, sampler)

particles.array[3,:] .= (xmax - xmin) ./ n_particles;

poisson = OneDPoisson( mesh )
kernel = ParticleMeshCoupling1D( mesh, n_particles, degree_smoother, :collocation)

ex = zeros(nx)
rho = zeros(nx)

for i_part = 1:particles.n_particles
    xi = particles.array[1, i_part]
    wi = particles.array[3, i_part]
    GEMPIC.add_charge!(rho, kernel, xi, wi)
end

compute_e_from_rho!(ex, poisson, rho)

problem = OneDPoissonPIC( poisson, kernel )

dt = 0.1
nsteps = 100
alpha = 0.1
kx = 0.5

propagator = OneDSplittingOperator( problem, particles ) 

energy = Float64[]

for j=1:nsteps

    ParticleInCell.operator_t!(propagator, 0.5dt)
    ParticleInCell.charge_deposition!(propagator)
    ParticleInCell.solve_fields!(propagator)
    ParticleInCell.operator_v!(propagator, dt)
    ParticleInCell.operator_t!(propagator, 0.5dt)

    push!(energy, compute_field_energy(problem))
          
end

t = collect(0:nsteps) .* dt
plot(log.(energy))

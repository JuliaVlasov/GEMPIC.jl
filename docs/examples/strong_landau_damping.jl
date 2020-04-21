# # Strong Landau Damping

# If we only use the subsystems $H_p$ and $H_E$, we are solving Vlasov--Poisson system.
# In this test, initial condition is as follows:
# ```math
# f_0(x,v) = \frac{1}{\sqrt{2\pi} \sigma}e^{-\frac{v^2}{2\sigma^2}}(1+\alpha \cos(kx)), \quad E_{10} (x) = \frac{\alpha}{k}\sin(kx), \ x \in [0,2\pi/k ),  \ v \in \mathbb{R}^2.
# ```

using GEMPIC, Plots

# The physical parameters are chosen as ``\sigma = 1``, ``k = 0.5``, 
# ``\alpha = 0.5``, and the numerical parameters as ``\Delta t = 0.05``, 
# ``n_x = 32`` and ``2 \times 10^5`` particles. We use second order Strang 
# splitting method in time. 

σ, μ = 1.0, 0.0
kx, α = 0.5, 0.5
xmin, xmax = 0, 2π/kx
∆t = 0.05
nx = 32 
n_particles = 200000
mesh = Mesh( xmin, xmax, nx)
spline_degree = 3

df = CosSumGaussian{1,2}([[kx]], [α], [[σ,σ]], [[μ,μ]] )

mass, charge = 1.0, 1.0
particle_group = ParticleGroup{1,2}(n_particles)   
sampler = ParticleSampler{1,2}( :sobol, true, n_particles)

sample!(particle_group, sampler, df, mesh)

xp = view(particle_group.array, 1, :)
vp = view(particle_group.array, 2:3, :)
wp = view(particle_group.array, 4, :);

p = plot(layout=(3,1))
histogram!(p[1,1], xp, weights=wp, normalize=true, bins = 100, lab = "")
plot!(p[1,1], x-> (1+α*cos(kx*x))/(2π/kx), 0., 2π/kx, lab="")
histogram!(p[2,1], vp[1,:], weights=wp, normalize=true, bins = 100, lab = "")
plot!(p[2,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
histogram!(p[3,1], vp[2,:], weights=wp, normalize=true, bins = 100, lab = "")
plot!(p[3,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")

# md savefig("histograms.svg"); nothing #hide
# md # ![](histograms.svg)

# Initialize the arrays for the spline coefficients of the fields

kernel_smoother1 = ParticleMeshCoupling( mesh, n_particles, spline_degree-1, :galerkin)    
kernel_smoother0 = ParticleMeshCoupling( mesh, n_particles, spline_degree, :galerkin)

# Initialize field solver



p = plot(layout=(1,2))
rho = zeros(Float64, nx)
efield_poisson = zeros(Float64, nx)
maxwell_solver = Maxwell1DFEM(mesh, spline_degree)
solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
xg = LinRange(xmin, xmax, nx)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, rho)
plot!(p[1,1], xg, sval, title=:rho, label=nothing)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, efield_poisson)
plot!(p[1,2], xg, sval, title=:ex, label=nothing )       
#md savefig("fields.svg"); nothing #hide
#md # ![](fields.svg)

# Simulation function

# You get better performance if your simulation is inside a function:

# +
function run( steps)

    σ, μ = 1.0, 0.0
    kx, α = 0.5, 0.5
    xmin, xmax = 0, 2π/kx
    dt = 0.01
    nx = 32 
    n_particles = 100000
    mesh = Mesh( xmin, xmax, nx)
    spline_degree = 3
    
    df = CosSumGaussian{1,2}([[kx]], [α], [[σ,σ]], [[μ,μ]] )
    
    mass, charge = 1.0, 1.0
    particle_group = ParticleGroup{1,2}(n_particles)   
    sampler = ParticleSampler{1,2}( :sobol, true, n_particles)
    
    sample!(particle_group, sampler, df, mesh)
    
    kernel_smoother1 = ParticleMeshCoupling( mesh, n_particles, spline_degree-1, :galerkin)    
    kernel_smoother0 = ParticleMeshCoupling( mesh, n_particles, spline_degree, :galerkin)
    
    rho = zeros(Float64, nx)
    efield_poisson = zeros(Float64, nx)

    maxwell_solver = Maxwell1DFEM(mesh, spline_degree)

    solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
    
    efield_dofs = [efield_poisson, zeros(Float64, nx)]
    bfield_dofs = zeros(Float64, nx)
        
    propagator = HamiltonianSplitting( maxwell_solver,
                                       kernel_smoother0, 
                                       kernel_smoother1, 
                                       particle_group,
                                       efield_dofs,
                                       bfield_dofs)
    
    efield_dofs_n = propagator.e_dofs
    
    thdiag = TimeHistoryDiagnostics( particle_group, maxwell_solver, 
                            kernel_smoother0, kernel_smoother1 )
    
    for j = 1:steps # loop over time
    
        strang_splitting!(propagator, dt, 1)
    
        solve_poisson!( efield_poisson, particle_group, 
                        kernel_smoother0, maxwell_solver, rho)
        
        write_step!(thdiag, j * dt, spline_degree, 
                        efield_dofs,  bfield_dofs,
                        efield_dofs_n, efield_poisson)
    end
    
    return thdiag.data

end
# -

@time results = run(5000) # change number of steps

plot(results[!,:Time], log.(results[!,:PotentialEnergyE1]))

# md savefig("potential_energy.svg"); nothing # hide
# md # ![](potential_energy.svg)

# *Time evolution of electric energy (semi-$\log_{10}$ scale).*

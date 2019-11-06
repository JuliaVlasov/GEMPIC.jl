
include("mesh.jl")
include("low_level_bsplines.jl")
include("splinepp.jl")
include("distributions.jl")
include("maxwell_1d_fem.jl")
include("particle_group.jl")
include("particle_mesh_coupling.jl")
include("particle_sampling.jl")
include("landau_damping.jl")
include("hamiltonian_splitting.jl")
include("hamiltonian_splitting_boris.jl")
include("diagnostics.jl")

const β = 0.0001
const k  = 1.25
const α  = 0.0
const σ  = [0.2,  0.005773502691896]
const μ  = [0.0, 0.0]

const nx   = 32
const xmin = 0.0
const xmax = 2π / k
const n_particles    = 100000
const sampling_case  = :sobol
const symmetric      = true
const splitting_case = :symplectic
const spline_degree  = 3

function run_simulation( steps )

    mesh   = Mesh( xmin, xmax, nx)
    domain = [xmin, xmax, xmax - xmin ]
    
    beta_sin_k(x) = beta * sin(2π * x / domain[3]) 
    
    df = CosSumGaussian{1,2}([[k]], [α], [σ], [μ] )
    
    v1min, v1max, nv1 = -1., 1., 64
    v2min, v2max, nv2 = -0.015, 0.015, 64
    v1 = LinRange(v1min, v1max, nv1) |> collect
    v2 = LinRange(v2min, v2max, nv2) |> collect
    f = zeros(Float64,(nv1,nv2))
    for i in eachindex(v1), j in eachindex(v2)
        f[i,j] = eval_v_density(df, [v1[i],v2[j]])
    end
    
     # Initialize the particles   (mass and charge set to 1.0 with one weight)
    mass, charge = 1.0, 1.0
    particle_group = ParticleGroup{1,2}( n_particles, mass, charge, 1)   
    
    sampler = ParticleSampler{1,2}( sampling_case, symmetric, n_particles)
    
    sample!(  particle_group, sampler, df, mesh)
    
    kernel_smoother1 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree-1, :galerkin)    
    kernel_smoother0 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree, :galerkin)
    rho = zeros(Float64, nx);
    efield_poisson = zeros(Float64, nx)
    # Initialize the field solver
    maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
    # efield by Poisson
    solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )    
    
    # Initialize the arrays for the spline coefficients of the fields
    efield_dofs = [copy(efield_poisson), zeros(Float64, nx)]
    bfield_dofs = zeros(Float64, nx)
        
    propagator = HamiltonianSplitting( maxwell_solver,
                                       kernel_smoother0, 
                                       kernel_smoother1, 
                                       particle_group,
                                       efield_dofs,
                                       bfield_dofs,
                                       domain);
    
    efield_dofs_n = propagator.e_dofs
    # -
    # bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))
    beta_cos_k(x) = β * cos(2π * x / domain[3]) 
    l2projection!( bfield_dofs, maxwell_solver, beta_cos_k, spline_degree-1)
    
    thdiag = TimeHistoryDiagnostics( particle_group, maxwell_solver, 
                            kernel_smoother0, kernel_smoother1 );
    
    Δt = 0.05
    
    for j = 1:steps # loop over time
    
        # Strang splitting
        strang_splitting!(propagator, Δt, 1)
    
        # Diagnostics
        solve_poisson!( efield_poisson, particle_group, 
                        kernel_smoother0, maxwell_solver, rho)
        
        write_step!(thdiag, j * Δt, spline_degree, 
                        efield_dofs,  bfield_dofs,
                        efield_dofs_n, efield_poisson)
    
    end

    return thdiag.data
    
end

run_simulation(1)
using Profile
Profile.clear()
@profile run_simulation(500)
Profile.print()

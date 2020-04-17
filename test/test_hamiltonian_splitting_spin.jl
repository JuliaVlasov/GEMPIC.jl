@testset "Hamiltonian splitting with spin" begin

    import GEMPIC: set_x, set_v
    import GEMPIC: set_weights, get_charge, add_charge!
    import GEMPIC: operatorHp1, operatorHp2, operatorHE, operatorHB

    # Tolerance for comparison of real numbers: set it here!
    eqv_tol = 1.0e-14

    # Set parameters
    n_particles     = 2
    eta_min         = 0.0
    eta_max         = 4π
    num_cells       = 10
    delta_t         = 0.1
    degree_smoother = 3
    rnd_seed        = 10

    mesh = Mesh( eta_min, eta_max, num_cells)

    pg = ParticleGroup{1,1}(n_particles, n_spin = 3)

    # Initialize kernel smoothers
    kernel_smoother_1 = ParticleMeshCoupling( mesh,
         n_particles, degree_smoother-1, :galerkin) 

    kernel_smoother_0 = ParticleMeshCoupling( mesh,
         n_particles, degree_smoother, :galerkin) 
    
    # Initialize Maxwell solver
    maxwell_solver = Maxwell1DFEM( mesh, degree_smoother)
    
    n_dofs = kernel_smoother_0.n_dofs

    ex         = zeros(Float64, n_dofs)
    efield     = ones( Float64, n_dofs)
    bfield     = ones( Float64, n_dofs)
    efield_ref = zeros(Float64, n_dofs)
    bfield_ref = zeros(Float64, n_dofs)
    rho        = zeros(Float64, n_dofs)

    for i_part = 1:n_particles
       xi = get_x(pg, i_part)[1]
       wi = get_charge( pg, i_part)
       add_charge!( rho, kernel_smoother_0, xi, wi)
    end

    @test true

end



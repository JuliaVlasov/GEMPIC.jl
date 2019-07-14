@testset "Hamiltonian splitting" begin

    import GEMPIC: set_common_weight, set_x, set_v
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

    domain = [eta_min, eta_max, eta_max - eta_min]

    pg = ParticleGroup{1,2}(n_particles, 1.0, 1.0, 1)

    set_common_weight( pg, 1.0)

    particle_info_ref = reshape([ 11.780972450961723, 
                                   5.4977871437821380,
                                  -1.5341205443525459, 
                                   0.15731068461017067,
                                   0.15731068461017067,
                                  -1.5341205443525459,
                                   6.8636759376074723,
                                   5.7026946767517002], n_particles, 4)

    # Initialize particles from particle_info_ref
    xi = zeros(2)
    for i_part = 1:n_particles
       xi[1]   = particle_info_ref[i_part, 1]
       set_x(pg, i_part, xi[1])
       xi[1:2] = particle_info_ref[i_part, 2:3]
       set_v(pg, i_part, xi)
       xi[1]   = particle_info_ref[i_part, 4]
       set_weights(pg, i_part, xi[1])
    end

    set_common_weight(pg, 1.0)

    # Initialize kernel smoothers
    kernel_smoother_1 = ParticleMeshCoupling( domain[1:2], [num_cells],
         n_particles, degree_smoother-1, :galerkin) 

    kernel_smoother_0 = ParticleMeshCoupling( domain[1:2], [num_cells],
         n_particles, degree_smoother, :galerkin) 
    
    # Initialize Maxwell solver
    maxwell_solver = Maxwell1DFEM( [eta_min, eta_max], num_cells,
                                    degree_smoother)
    
    n_dofs = kernel_smoother_0.n_dofs

    ex         = zeros(Float64, (n_dofs))
    efield_1   = ones( Float64, (n_dofs))
    efield_2   = ones( Float64, (n_dofs))
    bfield     = ones( Float64, (n_dofs))
    efield_ref = zeros(Float64, (n_dofs,2))
    bfield_ref = zeros(Float64, (n_dofs))
    rho        = zeros(Float64, (n_dofs))

    wi = zeros(1)
    for i_part = 1:n_particles
       xi = get_x(pg, i_part)
       wi = get_charge( pg, i_part)
       add_charge!( rho, kernel_smoother_0, xi, wi[1])
    end

    # Solve Poisson problem
    compute_e_from_rho!( efield_1, maxwell_solver, rho)

    propagator = HamiltonianSplitting( maxwell_solver,
         kernel_smoother_0, kernel_smoother_1, pg,
         [efield_1, efield_2], bfield, domain )

    operatorHp1(propagator, delta_t)

    # Compare to reference
    # Particle information after operatorV application 
    particle_info_ref = reshape([  11.627560396526469, 
                                    5.5135182122431550,
                                   -1.5341205443525459,
                                    0.15731068461017067,
                                    0.31072273904542569,
                                   -1.5498516128135633, 
                                    6.8636759376074723,
                                    5.7026946767517002], n_particles, 4)

    # Compare computed values to reference values
    for i_part=1:n_particles

       xi = get_x(pg, i_part)
 
       @test xi[1] ≈ particle_info_ref[i_part,1]

       xi = get_v(pg, i_part)
      
       @test xi[1] ≈ particle_info_ref[i_part,2]
       @test xi[2] ≈ particle_info_ref[i_part,3]

       xi = get_charge(pg, i_part)
       
       @test xi[1] ≈ particle_info_ref[i_part,4]

    end
    
    operatorHp2( propagator, delta_t)

    # Compare to reference
    # Particle information after operatorV application 
    particle_info_ref = reshape([ 11.627560396526469,
                                   5.5135182122431550,
                                  -1.5030482704480033,
                                  2.3255233288143329e-003,  
                                  0.31072273904542569,
                                  -1.5498516128135633, 
                                   6.8636759376074723,
                                   5.7026946767517002], n_particles, 4)

    # Compare computed values to reference values
    for i_part = 1:n_particles

       xi = get_x(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,1]

       xi = get_v(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,2]
       @test xi[2] ≈ particle_info_ref[i_part,3]

       xi[1] = get_charge(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,4]

    end

    operatorHE( propagator, delta_t)

    # Compare to reference
    # Particle information after operatorV application 
    particle_info_ref = reshape([ 11.627560396526469,
                                  5.5135182122431550,   
                                 -1.4311461960364338,
                                  2.7611337138493959e-002,  
                                  0.39511488788223620,
                                 -1.3903131848268702, 
                                  6.8636759376074723,
                                  5.7026946767517002 ], n_particles, 4)

    # Compare computed values to reference values
    for i_part = 1:n_particles

       xi = get_x(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,1]

       xi = get_v(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,2]
       @test xi[2] ≈ particle_info_ref[i_part,3]

       xi[1] = get_charge(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,4]

    end

    operatorHB(propagator, delta_t)

    # Compare to reference
    # Particle information after operatorV application 
    particle_info_ref = reshape([   11.627560396526469,
                                    5.5135182122431550,
                                   -1.4311461960364338,
                                    2.7611337138493959e-002,
                                    0.39511488788223620,
                                   -1.3903131848268702, 
                                    6.8636759376074723,
                                    5.7026946767517002], n_particles, 4)

    # Compare computed values to reference values
    for i_part=1:n_particles

       xi = get_x(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,1]

       xi = get_v(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,2]
       @test xi[2] ≈ particle_info_ref[i_part,3]

       xi = get_charge(pg, i_part)

       @test xi[1] ≈ particle_info_ref[i_part,4]

    end

    bfield_ref = [ 0.96241425319343565, 1.0704842857730226,
                   0.87749470807080665, 1.0585080035598995,
                   1.0312067899733146,  0.98524486372490983,
                   0.99163866889590213, 1.0359265946162064,
                   0.96576229888372844, 1.0213195333087741 ]

    efield_ref = reshape([ 0.32668967300827889,  
                           0.21111816665325256,
                          -3.2797353014133206, 
                           0.93182059869704992,
                           2.4099520529926513, 
                          -0.54074935301675420,
                          -6.6508143734355582e-002, 
                          -4.4526750013638647,
                           2.7314012140039408, 
                           2.4952249586987651,
                           1.2896538770270727,
                           0.45428428470806997,
                           1.9372031866314203,
                           1.2380946765497305,
                           0.83732795515227731, 
                           1.0244885213955153, 
                           1.1187987978241094,
                           0.68789307296286761,
                           1.0938576576751671,
                           0.85201507883794003], num_cells,2)

    @test maximum(abs.(bfield   .- bfield_ref))      ≈ 0.0 atol=1e-14
    @test maximum(abs.(efield_1 .- efield_ref[:,1])) ≈ 0.0 atol=1e-14
    @test maximum(abs.(efield_2 .- efield_ref[:,2])) ≈ 0.0 atol=1e-14

end



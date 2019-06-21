@testset "Hamiltonian splitting" begin

    import GEMPIC: set_common_weight, set_x, set_v
    import GEMPIC: set_weights, get_charge, add_charge

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

    pg = ParticleGroup{1,2}(n_particles, n_particles ,1.0, 1.0, 1)

    set_common_weight( pg, 1.0)

    particle_info_ref = reshape([11.780972450961723, 
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

    # Initialize kernel smoother    
    kernel_smoother_1 = ParticleMeshCoupling( domain[1:2], [num_cells],
         n_particles, degree_smoother-1, :galerkin) 
    kernel_smoother_0 = ParticleMeshCoupling( domain[1:2], [num_cells],
         n_particles, degree_smoother, :galerkin) 
    
    # Initialize Maxwell solver
    maxwell_solver = Maxwell1DFEM( [eta_min, eta_max], num_cells,
                                    degree_smoother)
    
    n_dofs = kernel_smoother_0.n_dofs
    ex         = zeros(Float64, (n_dofs))
    efield     = ones(Float64, (n_dofs,2))
    bfield     = ones(Float64, (n_dofs))
    efield_ref = zeros(Float64, (n_dofs,2))
    bfield_ref = zeros(Float64, (n_dofs))
    rho        = zeros(Float64, (n_dofs))
    rho_local  = zeros(Float64, (n_dofs))

    wi = zeros(1)
    for i_part = 1:n_particles
       xi = get_x(pg, i_part)
       wi[1] = get_charge( pg, i_part)
       add_charge( kernel_smoother_0, xi, wi[1], rho_local)
    end

    # Solve Poisson problem
    ex = zeros(Float64, n_dofs) 
    compute_e_from_rho( maxwell_solver, rho, ex)
    efield[:,1] .= ex

    propagator = HamiltonianSplitting( maxwell_solver,
         kernel_smoother_0, kernel_smoother_1, pg,
         efield, bfield, eta_min, eta_max-eta_min)

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

       xi[1] = get_charge(pg, i_part)
       
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

       xi[1:1] = get_charge(pg, i_part)

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

       xi[1:1] = get_charge(pg, i_part)

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

    @test maximum(abs.(bfield .- bfield_ref)) ≈ 0.0

    @test maximum(abs.(efield .- efield_ref)) ≈ 0.0

end


@testset "Hamiltonian splitting Boris" begin

    n_particles     = 2
    eta_min         = 0.0
    eta_max         = 4.0π
    num_cells       = 16
    delta_t         = 0.1_f64
    degree_smoother = 3

    pg = ParticleGroup{1,2}(n_particles, n_particles ,1.0, 1.0, 1)

    particle_info_ref = reshape( [11.780972450961723, 
                                 -1.5341205443525459,
                                  0.15731068461017067,
                                  6.8636759376074723,
                                  5.4977871437821380,
                                  0.15731068461017067,
                                 -1.5341205443525459,
                                  5.7026946767517002],4,n_particles)
    
    # Initialize particles from particle_info_ref
    xi = 0.0
    for i_part = 1:n_particles

       xi[1] = particle_info_ref[1, i_part]

       set_x(pg, i_part, xi)

       xi[1:2] = particle_info_ref[2:3, i_part]

       set_v(pg, i_part, xi)

       xi[1] = particle_info_ref(4, i_part)

       set_weights(pg, i_part, xi[1])

    end
    
    set_common_weight (pg, 1.0)

    # Initialize kernel smoother    

    kernel_smoother_1 = ParticleMeshCoupling( domain[1:2], [num_cells],
         n_particles, degree_smoother-1, :galerkin) 
    kernel_smoother_0 = ParticleMeshCoupling( domain(1:2), [num_cells],
         n_particles, degree_smoother, :galerkin) 
    
    maxwell_solver = MaxwellSolver( [eta_min, eta_max], num_cells, degree_smoother)
    
    efield      = zeros(Float64, (kernel_smoother_0.n_dofs,2))
    bfield      = zeros(Float64, (kernel_smoother_0.n_dofs))
    efield_ref  = zeros(Float64, (kernel_smoother_0.n_dofs,2))
    bfield_ref  = zeros(Float64, (kernel_smoother_0.n_dofs))
    rho         = zeros(Float64, (kernel_smoother_0.n_dofs))
    rho_local   = zeros(Float64, (kernel_smoother_0.n_dofs))

    efield = 1.0

    rho_local = 0.0
    for i_part = 1:n_particles

       xi = get_x( pg, i_part)

       wi[1] = get_charge(pg, i_part)

       add_charge(kernel_smoother_0, xi[1], wi[1], rho_local)

    end

    efield[:,1] .= compute_E_from_rho(rho)
    bfield      .= 1.0
    
    propagator = HamiltonSplittingBoris( maxwell_solver,
         kernel_smoother_0, kernel_smoother_1, pg,
         efield, bfield, eta_min, eta_max-eta_min)

    staggering( propagator, 0.5*delta_t )

    strang_splitting( propagator, delta_t, 1 )
    
    # Compare to reference
    # Particle information after operatorV application 
    particle_info_ref = reshape([ 11.590250917552364,
                                  -1.5236851980054489,
                                   0.40632305631363136,
                                   6.8636759376074723,
                                   5.5028810915315152,
                                   1.1611806341228641e-002,
                                  -1.4090781780562809,
                                   5.7026946767517002], 4,n_particles)

    for i_part=1:n_particles

         xi = get_x(pg, i_part)

         @test xi[1] ≈ particle_info_ref[1,i_part]

         xi = get_v(pg, i_part)

         @test xi[1] ≈ particle_info_ref[2,i_part]
         @test xi[2] ≈ particle_info_ref[3,i_part]

         xi[1:1] = get_charge(pg, i_part)

         @test xi[1] ≈ particle_info_ref[4,i_part]

    end

    bfield_ref = [1.0007152502482244,
                  0.99445857477120148,
                  1.0126078645236389,
                  0.97492340832308344,
                  1.0502493982082513,
                  0.88560153553873056,
                  1.1141170903091009,
                  0.95020604146500376, 
                  1.0248316646998006, 
                  0.98745124326634226,
                  1.0056520876592949,
                  0.99896954592982978,
                  0.99645963114519265,
                  1.0122347396384064, 
                  0.98742382899460557, 
                  1.0040980952792931 ]

    @test maximum(abs.(bfield .- bfield_ref)) ≈ 0.0 
    
    efield_ref = reshape( [ 1.7801632633330595,  
                            0.51717913695931361,
                            4.2412959490364575e-002, 
                           -1.1212332464258261, 
                           -1.1684637217416221,   
                           -3.7806966401897091,
                            3.7363585338991756,   
                            1.1932927327940017,
                            1.1033167281181586,    
                           -0.018180082909540971,
                           -0.56595990082616909,  
                           -1.6702327921823792,
                           -1.8045561507212267,
                           -4.3141455300350842,
                            4.8236042130256775,
                            1.5537560432220632,
                            0.98879315272141322,   
                            1.0323154046944396,
                            0.93329346828214133,    
                            1.1302445587553216,
                            0.73558670810944515,    
                            1.6340701469431378,
                            0.73779661553285170,    
                            1.1288774513495987,
                            0.93385001285635261,    
                            1.0324077177717137,
                            0.98801632510199811,   
                            0.99610949244376890,
                            1.0239154844069209,   
                            0.92782406399041817,
                            1.0265970800925086,   
                            0.99441071503466349 ], num_cells, 2)

    @test maximum(abs.(efield .- efield_ref))  ≈ 0.0
  
end 

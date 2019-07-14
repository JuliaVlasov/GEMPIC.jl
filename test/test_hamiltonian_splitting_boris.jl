import GEMPIC: set_x, set_v, set_weights, set_common_weight
import GEMPIC: get_charge, add_charge!
import GEMPIC: push_x_accumulate_j!

@testset "Hamiltonian splitting Boris" begin

    n_particles     = 2
    eta_min         = 0.0
    eta_max         = 4.0π
    num_cells       = 16
    delta_t         = 0.1
    degree_smoother = 3

    domain = [eta_min, eta_max, eta_max - eta_min]

    pg = ParticleGroup{1,2}(n_particles, 1.0, 1.0, 1)

    set_common_weight(pg, 1.0)

    particle_info_ref = reshape( [11.780972450961723, 
                                 -1.5341205443525459,
                                  0.15731068461017067,
                                  6.8636759376074723,
                                  5.4977871437821380,
                                  0.15731068461017067,
                                 -1.5341205443525459,
                                  5.7026946767517002],4,n_particles)
    
    # Initialize particles from particle_info_ref

    for i_part = 1:n_particles

       xi = particle_info_ref[1, i_part]
       set_x(pg, i_part, xi)
       xi = particle_info_ref[2:3, i_part]
       set_v(pg, i_part, xi)
       xi = particle_info_ref[4, i_part]
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
    
    efield  = [ ones(Float64, (kernel_smoother_0.n_dofs)),
                ones(Float64, (kernel_smoother_0.n_dofs))]

    bfield  = ones(Float64, (kernel_smoother_0.n_dofs))

    rho     = zeros(Float64, (kernel_smoother_0.n_dofs))

    for i_part = 1:n_particles

       xi = get_x( pg, i_part)

       wi = get_charge(pg, i_part)

       add_charge!( rho, kernel_smoother_0, xi, wi[1])

    end

    compute_e_from_rho!(efield[1], maxwell_solver, rho)

    propagator = HamiltonianSplittingBoris( maxwell_solver,
         kernel_smoother_0, kernel_smoother_1, pg,
         efield, bfield, domain)

    staggering!( propagator, 0.5*delta_t )

    strang_splitting!( propagator, delta_t, 1 )

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

         xi = get_charge(pg, i_part)

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

    @test maximum(abs.(efield[1] .- efield_ref[:,1])) ≈ 0.0 atol=1e-15
    @test maximum(abs.(efield[2] .- efield_ref[:,2])) ≈ 0.0 atol=1e-15
  
end

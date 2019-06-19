@testset "Hamiltonian splitting" begin

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

  pg = ParticleGroup(n_particles, n_particles ,1.0, 1.0, 1)

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
  xi = 0.0
  for i_part = 1:n_particles
     xi[1]   = particle_info_ref[i_part, 1]
     set_x(pg, i_part, xi)
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
  maxwell_solver = MaxwellSolver( [eta_min, eta_max], num_cells,
                                  degree_smoother)
  
  efield     = zeros(Float64, (kernel_smoother_0.n_dofs,2))
  bfield     = zeros(Float64, (kernel_smoother_0.n_dofs))
  efield_ref = zeros(Float64, (kernel_smoother_0.n_dofs,2))
  bfield_ref = zeros(Float64, (kernel_smoother_0.n_dofs))
  rho        = zeros(Float64, (kernel_smoother_0.n_dofs))
  rho_local  = zeros(Float64, (kernel_smoother_0.n_dofs))

  efield = 1.0

  rho_local = 0.0_f64
  for i_part = 1:n_particles
     xi = get_x(pg, i_part)
     wi[1] = get_charge( pg, i_part)
     add_charge( kernel_smoother_0, xi[1], wi[1], rho_local)
  end
  # Solve Poisson problem
  compute_E_from_rho( maxwell_solver, efield[:,1], rho)
  bfield = 1.0_f64

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
  end do
  
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
  end do

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

  n_particles     = 2
  eta_min         = 0.0
  eta_max         = 4.0π
  num_cells       = 16
  delta_t         = 0.1_f64
  degree_smoother = 3

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
     wi[1] = particle_group%get_charge( i_part)
     add_charge(kernel_smoother_0, xi[1], wi[1], rho_local)
  end

  call maxwell_solver%compute_E_from_rho(efield(:,1),&
       rho)
  bfield = 1.0_f64
  
  call propagator%init( maxwell_solver, &
       kernel_smoother_0, kernel_smoother_1, particle_group, &
       efield, bfield, &
       eta_min, eta_max-eta_min)
  call propagator%staggering( 0.5_f64*delta_t )

  call propagator%strang_splitting( delta_t, 1 )
  
  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([11.590250917552364_f64,       &
-1.5236851980054489_f64,       0.40632305631363136_f64,      &
  6.8636759376074723_f64,        5.5028810915315152_f64,     &
   1.1611806341228641E-002_f64,  -1.4090781780562809_f64,    &
    5.7026946767517002_f64 ], [4,n_particles])
  ! Compare computed values to reference values
    do i_part=1,n_particles
       xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(1,i_part))> EQV_TOL) then
        passed = .FALSE.
        print*, 'x(1) of particle ', i_part,' wrong'
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(2,i_part))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(1) of particle ', i_part,' wrong'
        print*, i_part, xi(1), particle_info_ref(i_part,2)
     elseif (abs(xi(2)-particle_info_ref(3,i_part))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(2) of particle ', i_part,' wrong'
        print*,   i_part, xi(2), particle_info_ref(3,i_part)
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(4,i_part))> EQV_TOL) then
        passed = .FALSE.
        print*, 'weight of particle ', i_part,' wrong'
     end if
  end do


  bfield_ref = [1.0007152502482244_f64,       0.99445857477120148_f64,   &
     1.0126078645236389_f64,       0.97492340832308344_f64,        &
1.0502493982082513_f64,       0.88560153553873056_f64,        &
1.1141170903091009_f64,       0.95020604146500376_f64,        &
1.0248316646998006_f64,       0.98745124326634226_f64,        &
1.0056520876592949_f64,       0.99896954592982978_f64,       &
0.99645963114519265_f64,        1.0122347396384064_f64,       &
0.98742382899460557_f64,        1.0040980952792931_f64   ]
  error = maxval(abs(bfield-bfield_ref))  
 if (error > eqv_tol ) then
     print*, 'bfield error too large.'
     passed = .false.
  end if
  
  efield_ref = reshape( [1.7801632633330595_f64,  0.51717913695931361_f64, &
       4.2412959490364575E-002_f64,  -1.1212332464258261_f64,    &
       -1.1684637217416221_f64,       -3.7806966401897091_f64,&
       3.7363585338991756_f64,        1.1932927327940017_f64,&
       1.1033167281181586_f64,       -0.018180082909540971_f64,&
       -0.56595990082616909_f64,       -1.6702327921823792_f64, &
       -1.8045561507212267_f64,       -4.3141455300350842_f64,  &
       4.8236042130256775_f64,        1.5537560432220632_f64,   &
       0.98879315272141322_f64,        1.0323154046944396_f64,   &
       0.93329346828214133_f64,        1.1302445587553216_f64,  &
       0.73558670810944515_f64,        1.6340701469431378_f64,  &
       0.73779661553285170_f64,        1.1288774513495987_f64,   &
       0.93385001285635261_f64,       1.0324077177717137_f64,    &
       0.98801632510199811_f64,       0.99610949244376890_f64,  &
       1.0239154844069209_f64,       0.92782406399041817_f64,    &
       1.0265970800925086_f64,       0.99441071503466349_f64   ],&
       [num_cells, 2])

  error = maxval(abs(efield-efield_ref))  
  if (error > eqv_tol ) then
     print*, 'efield error too large.'
     passed = .false.
  end if
  
end 

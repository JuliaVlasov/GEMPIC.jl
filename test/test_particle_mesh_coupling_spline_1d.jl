using GEMPIC

import GEMPIC: set_common_weight, set_x, set_v, set_weights

function test_particle_mesh_coupling_spline_1d()

    n_cells       = 10            # Number of cells
    n_particles   = 4             # Number of particles
    spline_degree = 3             # Spline degree
    domain        = [0.0, 2.0]    # x_min, x_max
    x_vec = [0.1, 0.65, 0.7, 1.5] # Particle positions
    v_vec = [1.5  -3.0  0.0  6.0; 
             0.0   0.5  0.0  0.0]'

    particle_group = ParticleGroup1D2V( n_particles, n_particles ,1.0, 1.0, 1)
  
    set_common_weight(particle_group, 1.0/n_particles)

    for i_part = 1:n_particles
        set_x(particle_group, i_part, x_vec[i_part])
        set_weights(particle_group, i_part, 1.0)
        set_v(particle_group, i_part, v_vec[i_part,:])
    end
  
    values_grid = zeros(Float64,(4,1,4))

    values_grid[:,1,1] .= [ 2.0833333333333332e-002,  
                            0.47916666666666663, 
                            0.47916666666666663, 
                            2.0833333333333332e-002]

    values_grid[:,1,3] .= values_grid[:,1,1]   

    values_grid[:,1,4] .= values_grid[:,1,1] 

    values_grid[:,1,2] .= [7.0312500000000000e-002,  
                           0.61197916666666663,
                           0.31510416666666663,        
                           2.6041666666666665E-003 ]

    kernel = ParticleMeshCoupling( domain, [n_cells], n_particles, 
                 spline_degree, :collocation)
  
#=
  ! Check that the constructors for the abstract type are working.
  call sll_s_new_particle_mesh_coupling_spline_1d_ptr(ksp, domain, [n_cells], n_particles, spline_degree, sll_p_collocation)

  ! Accumulate rho
  rho_dofs = 0.0
  rho_dofs1 = 0.0
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call kernel%add_charge(xi(1), wi(1), rho_dofs)
     call kernel%add_charge_pp(xi(1), wi(1), rho_dofs1)
  end do
  !rho_dofs = rho_dofs
  rho_dofs_ref = 0.0
  rho_dofs_ref(8:10) = values_grid(1:3,1,1)
  rho_dofs_ref(1) = values_grid(4,1,1)
  rho_dofs_ref(1:4) = rho_dofs_ref(1:4) + values_grid(:,1,2) + values_grid(:,1,3)
  rho_dofs_ref(5:8) = rho_dofs_ref(5:8) + values_grid(:,1,4)
  rho_dofs_ref = rho_dofs_ref/real(n_particles, f64) * real(n_cells,f64)/domain(2)
  error = maxval(abs(rho_dofs-rho_dofs_ref))
  error1 = maxval(abs(rho_dofs1-rho_dofs_ref))

  if (error > 1.e-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge .'
  elseif (error1 >1.e-14)then 
     passed = .FALSE.
     print*, 'Error in procedure add_charge_pp .'
  end if


  ! Test j accumulations
  j_dofs = 0.0
  j_dofs1 = 0.0
  b_dofs = 0.0
  do i_part=1, n_particles 
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     vi = particle_group%get_v(i_part)
     vi1 = vi
     x_new = xi(1) + vi(1)/10.0

     call kernel%add_current_update_v( xi(1), x_new, wi(1), 1.0, b_dofs, vi, j_dofs )
     call kernel%add_current_update_v_pp( xi(1), x_new, wi(1), 1.0, b_dofs, vi1, j_dofs1 )
     
  end do
  j_dofs_ref = [ 2.4617513020833336D-002,   4.0690104166666692D-005, 0.0, &
       0.0, 0.0, 0.0,  0.0,  6.5104166666666674D-004, &
       4.6834309895833329D-002,    0.11535644531249999];
  j_dofs_ref = j_dofs_ref + [ -0.16219075520833331, -0.16219075520833331, &
       -2.52685546875D-002,  -4.0690104166666692D-005 , &
       0.0, 0.0, 0.0, 0.0,  &
       -4.0690104166666692D-005,  -2.52685546875D-002];
  j_dofs_ref = j_dofs_ref + [ 6.5104166666666696D-004,  0.0, 0.0,  &
       0.0, 6.5104166666666674D-004, 5.0130208333333329D-002,  &
       0.19986979166666666, 0.24869791666666663,  &
       0.19986979166666666, 5.0130208333333329D-002];
 
  error = maxval(abs(j_dofs-j_dofs_ref))
  error1= maxval(abs(j_dofs1-j_dofs_ref))
  
  if (error > 1.e-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_current_update_v.'

  elseif (error1 > 1.e-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_current_update_v_spline_pp.'
  end if


  call sll_s_spline_pp_b_to_pp_1d(kernel%spline_pp,n_cells,rho_dofs,rho_dofs_pp)
   
  ! Test function evaluation
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     call kernel%evaluate(xi(1), rho_dofs, particle_values(i_part))
     call kernel%evaluate_pp(xi(1), rho_dofs_pp, particle_values1(i_part))
  end do
  particle_values_ref = [1.1560058593749998,       2.3149278428819446, &
       2.2656250000000000,        1.1512586805555554]/domain(2);
  error = maxval(abs(particle_values-particle_values_ref))
  error1 = maxval(abs(particle_values1-particle_values_ref))
  !print*,'fehler=', maxval(abs(particle_values1-particle_values))

  if (error > 1.e-14) then
     passed = .FALSE.
     print*, 'Error in procedure evaluate_field_single.'
  else  if (error1 > 1.e-14) then
     passed = .FALSE.
     print*, 'Error in procedure evaluate_field_single_spline_pp.'
  end if

  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  call kernel%free()
  call particle_group%free()


=#
end

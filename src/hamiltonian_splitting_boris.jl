"""
Boris pusher in GEMPIC framework (spline finite elements)
Reference: Kraus, Kormann, Sonnendrücker, Morrison: GEMPIC: Geometric ElectroMagnetic Particle-In-Cell Methods

Solves Vlasov-Maxwell with PIC and spline finite elements with Boris pusher

- DoFs describing the magnetic field at time t_{n+1/2} (used for push)
- DoFs for kernel representation of current density. 
"""
struct HamiltonianSplittingBoris <: AbstractSplitting

     maxwell_solver    :: AbstractMaxwellSolver
     kernel_smoother_0 :: KernelSmoother #(order p+1)
     kernel_smoother_1 :: KernelSmoother #(order p)
     particle_group    :: ParticleGroup
     spline_degree     :: Int64
     Lx                :: Float64
     x_min             :: Float64 # Lower bound for x domain
     delta_x           :; Float64 # Grid spacing

     cell_integrals_0  :: SVector
     cell_integrals_1  :: SVector

     efield_dofs       :: Array{Float64, 2}
     efield_dofs_mid   :: Array{Float64, 2}
     bfield_dofs       :: Array{Float64, 2}
     bfield_dofs_mid   :; Array{Float64, 1}
     j_dofs            :: Array{Float64, 2}

end 

contains

"""
Second order Boris pusher using staggered grid
- self : time splitting object 
- dt   : time step
- number_steps : number of time steps
"""
function operator_boris(self, dt, number_steps)

    for i_step = 1:number_steps
        # (1) Compute B_{n+1/2} from B_n
        bfield_dofs_mid = bfield_dofs
        compute_B_from_E( maxwell_solver,
                          dt, 
                          self.efield_dofs_mid[:,2], 
                          self.bfield_dofs)

       self%bfield_dofs_mid = (self%bfield_dofs_mid + self%bfield_dofs)*0.5_f64
       
       ! (2) Propagate v: v_{n-1/2} -> v_{n+1/2}
       ! (2a) Half time step with E-part
       call push_v_epart( self, dt*0.5_f64 )
       ! (2b) Full time step with B-part
       call push_v_bpart( self, dt )
       ! (2c) Half time step with E-part
       call push_v_epart( self, dt*0.5_f64 )
       
       ! (3) Propagate x: x_n -> x_{n+1}. Includes also accumulation of j_x, j_y
       call push_x_accumulate_j ( self, dt )
       
       ! (4) Compute E_{n+1}       
       self%efield_dofs = self%efield_dofs_mid
       self%j_dofs = dt*self%j_dofs
       ! Ex part:
       call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, self%efield_dofs_mid(:,1))
       ! TODO: Combine the two steps for efficiency
       ! Ey part:    
       call self%maxwell_solver%compute_E_from_B(&
            dt, self%bfield_dofs, self%efield_dofs_mid(:,2))
       call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs_mid(:,2))
       
    end do

  end subroutine operator_boris

  !> Pusher for E \nabla_v part
  subroutine push_v_epart (self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v_boris), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_new(3), xi(3)
    sll_real64 :: efield(2)
    sll_real64 :: qm


    qm = self%particle_group%species%q_over_m();
    ! V_new = V_old + dt * E
    do i_part=1,self%particle_group%n_particles
       ! Evaluate efields at particle position
       xi = self%particle_group%get_x(i_part)
       call self%kernel_smoother_1%evaluate &
            (xi(1), self%efield_dofs_mid(:,1), efield(1))
       call self%kernel_smoother_0%evaluate &
            (xi(1), self%efield_dofs_mid(:,2), efield(2))
       v_new = self%particle_group%get_v(i_part)
       v_new(1:2) = v_new(1:2) + dt* qm * efield
       call self%particle_group%set_v(i_part, v_new)
    end do
    

  end subroutine push_v_Epart

  !>  Pusher for vxB part
   subroutine push_v_bpart (self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v_boris), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: vi(3), v_new(3), xi(3)
    sll_real64 :: bfield, M11, M12
    sll_real64 :: qmdt

    qmdt = self%particle_group%species%q_over_m()*0.5_f64*dt;
    
    do i_part=1,self%particle_group%n_particles
       vi= self%particle_group%get_v(i_part)
       xi = self%particle_group%get_x(i_part)
       call self%kernel_smoother_1%evaluate &
            (xi(1), self%bfield_dofs_mid, bfield)

       bfield = qmdt*bfield
       M11 = 1.0_f64/(1.0_f64 + bfield**2) 
       M12 = M11*bfield*2.0_f64
       M11 = M11*(1-bfield**2)

       v_new(1) = M11 * vi(1) + M12 * vi(2)
       v_new(2) = - M12 * vi(1) + M11 * vi(2)
       v_new(3) = 0.0_f64
       
       call self%particle_group%set_v(i_part, v_new)

    end do

  end subroutine push_v_bpart

  !> Pusher for x and accumulate current densities
  subroutine push_x_accumulate_j (self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v_boris), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: x_new(3), vi(3), wi(1), x_old(3)
    sll_int32  :: n_cells
    sll_real64 :: qoverm


    n_cells = self%kernel_smoother_0%n_dofs


    self%j_dofs_local = 0.0_f64

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  X_new = X_old + dt * V
    do i_part=1,self%particle_group%n_particles  
       ! Read out particle position and velocity
       x_old = self%particle_group%get_x(i_part)
       vi = self%particle_group%get_v(i_part)

       ! Then update particle position:  X_new = X_old + dt * V
       x_new = x_old + dt * vi

       ! Get charge for accumulation of j
       wi = self%particle_group%get_charge(i_part)
       qoverm = self%particle_group%species%q_over_m();

       ! TODO: Check here also second posibility with sum of two accumulations
       ! Accumulate jx
       call self%kernel_smoother_1%add_charge( [(x_old(1)+x_new(1))*0.5_f64], wi(1)*vi(1), &
            self%j_dofs_local(:,1))
       ! Accumulate jy
       call self%kernel_smoother_0%add_charge( [(x_old(1)+x_new(1))*0.5_f64], wi(1)*vi(2), &
            self%j_dofs_local(:,2))
      
       x_new(1) = modulo(x_new(1), self%Lx)
       call self%particle_group%set_x(i_part, x_new)

    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j_dofs(:,1))
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
         n_cells, MPI_SUM, self%j_dofs(:,2))
    

  end subroutine push_x_accumulate_j

 !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_1d2v_boris(&
       self, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx) 
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v_boris), intent(out) :: self !< time splitting object 
    class(sll_c_maxwell_1d_base), pointer,          intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling), pointer,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling), pointer,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_c_particle_group_base), pointer,      intent(in)  :: particle_group !< Particle group
    sll_real64, pointer,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, pointer,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.

    !local variables
    sll_int32 :: ierr

    self%maxwell_solver => maxwell_solver
    self%kernel_smoother_0 => kernel_smoother_0
    self%kernel_smoother_1 => kernel_smoother_1
    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs

#ifndef __PGI
    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_dofs == self%kernel_smoother_1%n_dofs )
#endif
    SLL_ALLOCATE(self%j_dofs(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%bfield_dofs_mid(self%kernel_smoother_1%n_dofs), ierr)
    SLL_ALLOCATE(self%efield_dofs_mid(self%kernel_smoother_1%n_dofs,2), ierr)

    self%spline_degree = 3
    self%x_min = x_min
    self%Lx = Lx
    self%delta_x = self%Lx/self%kernel_smoother_1%n_dofs
    
    self%cell_integrals_1 = [0.5_f64, 2.0_f64, 0.5_f64]
    self%cell_integrals_1 = self%cell_integrals_1 / 3.0_f64

    self%cell_integrals_0 = [1.0_f64,11.0_f64,11.0_f64,1.0_f64]
    self%cell_integrals_0 = self%cell_integrals_0 / 24.0_f64


  end subroutine initialize_pic_vm_1d2v_boris

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v_boris(self)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v_boris), intent( inout ) :: self !< time splitting object 

    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    deallocate(self%bfield_dofs_mid)
    deallocate(self%efield_dofs_mid)
    self%maxwell_solver => null()
    self%kernel_smoother_0 => null()
    self%kernel_smoother_1 => null()
    self%particle_group => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

  end subroutine delete_pic_vm_1d2v_boris


  !---------------------------------------------------------------------------!
  !> Propagate E_0 to E_{1/2} and x_0 to x_{1/2} to initialize the staggering
  subroutine staggering_pic_vm_1d2v_boris(self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v_boris), intent( inout ) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    call push_x_accumulate_j (self, dt*0.5_f64)

    ! (4) Compute E_{n+1}       
    self%efield_dofs_mid = self%efield_dofs
    self%j_dofs = dt*0.5_f64*self%j_dofs
    ! Ex part:
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, self%efield_dofs_mid(:,1))
    ! TODO: Combine the two steps for efficiency
    ! Ey part:    
    call self%maxwell_solver%compute_E_from_B(&
         dt*0.5_f64, self%bfield_dofs, self%efield_dofs_mid(:,2))
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs_mid(:,2))

  end subroutine staggering_pic_vm_1d2v_boris


end module sll_m_hamiltonian_splitting_pic_vm_1d2v_boris

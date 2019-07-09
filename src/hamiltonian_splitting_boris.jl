"""
Boris pusher in GEMPIC framework (spline finite elements)
Reference: Kraus, Kormann, Sonnendrücker, Morrison: GEMPIC: Geometric ElectroMagnetic Particle-In-Cell Methods

Solves Vlasov-Maxwell with PIC and spline finite elements with Boris pusher

- DoFs describing the magnetic field at time t_{n+1/2} (used for push)
- DoFs for kernel representation of current density. 
- `maxwell_solver`    : Maxwell solver
- `kernel_smoother_0` : Kernel smoother
- `kernel_smoother_1` : Kernel smoother
- `particle_group`    : Particle group
- `efield_dofs`       : array for the coefficients of the efields 
- `bfield_dofs`       : array for the coefficients of the bfield
- `x_min`             : Lower bound of x domain
- `Lx`                : Length of the domain in x direction.

"""
struct HamiltonianSplittingBoris <: AbstractSplitting

     maxwell_solver    :: AbstractMaxwellSolver
     kernel_smoother_0 :: KernelSmoother #(order p+1)
     kernel_smoother_1 :: KernelSmoother #(order p)
     particle_group    :: ParticleGroup
     spline_degree     :: Int64
     Lx                :: Float64
     x_min             :: Float64 # Lower bound for x domain
     delta_x           :: Float64 # Grid spacing

     cell_integrals_0  :: SVector
     cell_integrals_1  :: SVector

     e_dofs_1     :: Array{Float64, 1}
     e_dofs_2     :: Array{Float64, 1}
     e_dofs_1_mid :: Array{Float64, 1}
     e_dofs_2_mid :: Array{Float64, 1}
     b_dofs       :: Array{Float64, 1}
     b_dofs_mid   :: Array{Float64, 1}
     j_dofs_1     :: Array{Float64, 1}
     j_dofs_2     :: Array{Float64, 1}

     function HamiltonianSplittingBoris( maxwell_solver,
                                         kernel_smoother_0,
                                         kernel_smoother_1,
                                         particle_group,
                                         efield_dofs,
                                         bfield_dofs,
                                         x_min,
                                         Lx ) 

         e_dofs_1_mid = zeros(Float64, kernel_smoother_1.n_dofs)
         e_dofs_2_mid = zeros(Float64, kernel_smoother_1.n_dofs)
         j_dofs_1     = zeros(Float64, kernel_smoother_0.n_dofs)
         j_dofs_2     = zeros(Float64, kernel_smoother_0.n_dofs)
         b_dofs_mid   = zeros(Float64, kernel_smoother_1.n_dofs)

         spline_degree = 3
         delta_x       = Lx/kernel_smoother_1.n_dofs
    
         cell_integrals_1 = SVector([0.5,  2.0,  0.5] ./ 3.0)
         cell_integrals_0 = SVector([1.0, 11.0, 11.0, 1.0] ./ 24.0)

         new( maxwell_solver, kernel_smoother_0, kernel_smoother_1,
              particle_group, spline_degree, Lx, x_min, delta_x,
              cell_integrals_0, cell_integrals_1,
              e_dofs_1, e_dofs_2, e_dofs_1_mid, e_dofs_2_mid,
              b_dofs, b_dofs_mid, j_dofs_1, j_dofs_2 )

    end

end

"""
Propagate ``E_0`` to ``E_{1/2}`` and ``x_0`` to ``x_{1/2}`` to initialize 
the staggering
- self : time splitting object 
- dt   : time step
"""
function staggering_pic_vm_1d2v_boris(self, dt)

    push_x_accumulate_j!(self, dt*0.5)

    # (4) Compute E_{n+1}
    self.e_dofs_1_mid .= self.e_dofs_1
    self.e_dofs_2_mid .= self.e_dofs_2
    self%j_dofs = dt*0.5_f64*self%j_dofs
    # Ex part:
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, self%efield_dofs_mid(:,1))
    ! TODO: Combine the two steps for efficiency
    ! Ey part:    
    call self%maxwell_solver%compute_E_from_B(&
         dt*0.5_f64, self%bfield_dofs, self%efield_dofs_mid(:,2))
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs_mid(:,2))

end subroutine staggering_pic_vm_1d2v_boris


"""
    operator_boris(self, dt, number_steps)

Second order Boris pusher using staggered grid
- self : time splitting object 
- dt   : time step
- number_steps : number of time steps
"""
function operator_boris(self, dt, number_steps)

    for i_step = 1:number_steps

        # (1) Compute B_{n+1/2} from B_n
        bfield_dofs_mid = bfield_dofs
        compute_b_from_e!(self.b_dofs,
                          maxwell_solver,
                          dt, 
                          self.e_dofs_2_mid) 

        self.b_dofs_mid .+= self.b_dofs*0.5
       
        # (2) Propagate v: v_{n-1/2} -> v_{n+1/2}
        # (2a) Half time step with E-part
        push_v_epart!( self, dt*0.5 )
        # (2b) Full time step with B-part
        push_v_bpart!( self, dt )
        # (2c) Half time step with E-part
        push_v_epart!( self, dt*0.5 )
        
        # (3) Propagate x: x_n -> x_{n+1}. Includes also accumulation of j_x, j_y
        push_x_accumulate_j!( self, dt )
        
        # (4) Compute E_{n+1}       
        self.e_dofs_1 .= self.e_dofs_1_mid
        self.e_dofs_2 .= self.e_dofs_2_mid
        self.j_dofs_1 .= dt * self.j_dofs_1
        self.j_dofs_2 .= dt * self.j_dofs_2

        # Ex part:
        compute_e_from_j!(self.e_dofs_1_mid, maxwell_solver, self.j_dofs_1, 1)

        # Ey part:    
        compute_e_from_b!(self.e_dofs_2_mid, maxwell_solver, dt, self.b_dofs) 
        compute_e_from_j!(self.e_dofs_2_mid, maxwell_solver, self.j_dofs_2_dofs, 2)
       
    end

end

"""
Pusher for ``E \\nabla_v `` part
"""
function push_v_epart(self, dt)

    qm = self.particle_group.q_over_m

    # V_new = V_old + dt * E
    for i_part = 1:self.particle_group.n_particles

        # Evaluate efields at particle position
        xi = self.particle_group.get_x(i_part)

        efield[1] = evaluate(self.kernel_smoother_1, xi[1], self.e_dofs_1_mid)
        efield[2] = evaluate(self.kernel_smoother_0, xi[1], self.e_dofs_2_mid)

        v_new  = get_v(self.particle_group, i_part)
        v_new .= v_new .+ dt * qm * efield

        set_v!(particle_group, i_part, v_new)

    end
    

end

"""
  Pusher for vxB part
"""
function push_v_bpart(self, dt)

    qmdt = self.particle_group.q_over_m * 0.5 * dt
    
    do i_part=1:self.particle_group.n_particles

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

    end

end

"""
Pusher for x and accumulate current densities
"""
function push_x_accumulate_j (self, dt)

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
    

end 


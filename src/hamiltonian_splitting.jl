using StaticArrays

"""
Hamiltonian splitting type for Vlasov-Maxwell
- Integral over the spline function on each interval (order p+1)
- Integral over the spline function on each interval (order p)
!< DoFs describing the two components of the electric field
DoFs describing the magnetic field
DoFs for kernel representation of current density. 
MPI-processor local part of one component of \a j_dofs
"""

struct HamiltonianSplitting

    maxwell_solver    :: AbstractMaxwellSolver
    kernel_smoother_0 :: ParticleMeshCoupling
    kernel_smoother_1 :: ParticleMeshCoupling
    particle_group    :: ParticleGroup

    spline_degree     :: Int64
    Lx                :: Float64
    x_min             :: Float64
    delta_x           :: Float64

    cell_integrals_0  :: SVector
    cell_integrals_1  :: SVector

    efield_dofs       :: Array{Float64,2}
    bfield_dofs       :: Array{Float64,2}
    j_dofs            :: Array{Float64,2}
    j_dofs_local      :: Array{Float64,2}

    function HamiltonianSplitting( maxwell_solver,
                                   kernel_smoother_0,
                                   kernel_smoother_1,
                                   particle_group,
                                   efield_dofs,
                                   bfield_dofs,
                                   x_min,
                                   Lx) 

        # Check that n_dofs is the same for both kernel smoothers.
        @assert kernel_smoother_0.n_dofs == kernel_smoother_1.n_dofs

        j_dofs       = zeros(Float64,(kernel_smoother_0.n_dofs,2))
        j_dofs_local = zeros(Float64,(self.kernel_smoother_0.n_dofs,2))

        spline_degree = 3
        delta_x = Lx/kernel_smoother_1.n_dofs
        
        cell_integrals_1 = SVector{3}([0.5, 2.0, 0.5] ./ 3.0)
        cell_integrals_0 = SVector{4}([1.0,11.0,11.0,1.0] ./ 24.0)

        new( maxwell_solver, kernel_smoother_0, kernel_smoother_1, 
             particle_group, spline_degree, Lx, x_min, delta_x,
             cell_integrals_0, cell_integrals_1,
             efield_dofs, bfield_dofs, j_dofs, j_dofs_local)

    end

end

"""
Strang splitting
- time splitting object 
- time step
- number of time steps
"""
function strang_splitting( self         :: HamiltonianSplitting,
                           dt           :: Float64, 
                           number_steps :: Int64)

    for i_step = 1:number_steps
       operatorHB(  0.5dt)
       operatorHE(  0.5dt)
       operatorHp2( 0.5dt)
       operatorHp1( 1.0dt)
       operatorHp2( 0.5dt)
       operatorHE(  0.5dt)
       operatorHB(  0.5dt)
    end

end 

"""
Lie splitting
"""
function lie_splitting(self         :: HamiltonianSplitting,
                       dt           :: Float64, 
                       number_steps :: Int64)

    for i_step = 1:number_steps
       operatorHE(dt)
       operatorHB(dt)
       operatorHp1(dt)
       operatorHp2(dt)
    end

end 

"""
Lie splitting (oposite ordering)
"""
function lie_splitting_back(self         :: HamiltonianSplitting,
                            dt           :: Float64, 
                            number_steps :: Int64)

    for i_step = 1:number_steps
       operatorHp2(dt)
       operatorHp1(dt)
       operatorHB(dt)
       operatorHE(dt)
    end

end
  

"""
Push Hp1: Equations to solve are
```math
\\begin{eqnaray}
\\partial_t f + v_1 \\partial_{x_1} f = 0 & -> &  X_new = X_old + dt V_1 \\\\
V_new,2 = V_old,2 + \\int_0 h V_old,1 B_old
\\partial_t E_1 = - \\int v_1 f(t,x_1, v) dv & -> & E_{1,new} = E_{1,old} - 
\\int \\int v_1 f(t,x_1+s v_1,v) dv ds \\partial_t E_2 = 0 & -> & E_{2,new} = E_{2,old} \\\
\\partial_t B = 0 & => & B_new = B_old 
\\end{eqnarray}

 Here we have to accumulate j and integrate over the time interval.
 At each k=1,...,n_grid, we have for s \in [0,dt]:
 j_k(s) =  \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h) v_k,
 where h is the grid spacing and N the normalized B-spline
 In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
 j_k(s) = \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i) v_k.
 Now, we want the integral 
 \int_{0..dt} j_k(s) d s = \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i) ds =  \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt}  N(y_k + w v_{1,k}-y_i) dw

For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.

Then update particle position:  X_new = X_old + dt * V
"""
function operatorHp1(self :: HamiltonianSplitting, dt :: Float64)

    self :: HamiltonianSplitting 
    dt   :: Float64 

    n_cells = kernel_smoother_0.n_dofs

    j_dofs_local = 0.0

    for i_part=1:self.particle_group.n_particles  
       # Read out particle position and velocity
       x_old = self.particle_group.get_x(i_part)
       vi = self%particle_group%get_v(i_part)

       # Then update particle position:  X_new = X_old + dt * V
       x_new = x_old + dt * vi

       # Get charge for accumulation of j
       wi = self.particle_group.get_charge(i_part)
       qoverm = self.particle_group.species%q_over_m();

       kernel_smoother_1.add_current_update_v( x_old, x_new, wi[1],
            qoverm, self.bfield_dofs, vi, self.j_dofs_local[:,1])
       # Accumulate rho for Poisson diagnostics
       self.kernel_smoother_0.add_charge( x_new, wi[1], self.j_dofs_local[:,2])
      
       x_new[1] = modulo(x_new[1], self.Lx)
       particle_group.set_x(i_part, x_new)
       particle_group.set_v(i_part, vi)

    end

    self.j_dofs = 0.0

    # Update the electric field.
    compute_E_from_j(maxwell_solver, self.j_dofs[:,1], 1, self.efield_dofs[:,1])

end

"""
Push Hp2: Equations to solve are

```math
\\begin{eqnarray}
X_new  =  X_old && \\\\
V_new,1 = V_old,1 + \\int_0 h V_old,2 B_old && \\\\
\\partial_t E_1 = 0 & -> & E_{1,new} = E_{1,old}  \\\\
\\partial_t E_2 = - \\int v_2 f(t,x_1, v) dv & -> & \\\\ E_{2,new} = E_{2,old} 
- \\int \\int v_2 f(t,x_1+s v_1,v) dv ds
\\partial_t B = 0 => B_new = B_old
```
"""
function operatorHp2(self, dt)
    
    n_cells = self.kernel_smoother_0.n_dofs

    j_dofs_local = 0.0

    qm = self.particle_group.species.q_over_m()
    # Update v_1
    for i_part=1:self.particle_group.n_particles

       # Evaluate bfield at particle position (splines of order p)
       xi = particle_group.get_x(i_part)
       evaluate(kernel_smoother_1, xi[1], self.bfield_dofs, bfield)
       vi = particle_group.get_v(i_part)
       vi[1] = vi[1] + dt*qm*vi[2]*bfield
       set_v(self.particle_group, i_part, vi)

       xi = particle_group.get_x(i_part)

       # Scale vi by weight to combine both factors 
       #for accumulation of integral over j

       wi = get_charge(self.particle_group, i_part)*vi[2]

       add_charge(kernel_smoother_0, xi[1:1], wi[1], 
                  self.j_dofs_local[:,2]) 

    end

    self.j_dofs = 0.0

    # Update the electric field. Also, we still need to scale with 1/Lx 
    
    self.j_dofs[:,2] = self.j_dofs[:,2] * dt

    compute_E_from_j( maxwell_solver, 
                      self.j_dofs[:,2], 2, 
                      self.efield_dofs[:,2])
    
end
  
"""
Push H_E: Equations to be solved
```math
\\begin{eqnarray}
\\partial_t f + E_1 \\partial_{v_1} f + E_2 \\partial_{v_2} f = 0 &->& V_new = V_old + dt * E \\\\
\\partial_t E_1 = 0 &->& E_{1,new} = E_{1,old} \\\\
\\partial_t E_2 = 0 &->& E_{2,new} = E_{2,old} \\\\
\\partial_t B + \\partial_{x_1} E_2 = 0 &->& B_new = B_old - dt \\partial_{x_1} E_2
\\end{eqnarray}
```
"""
function operatorHE(self :: HamiltonianSplitting, dt)

    qm = self.particle_group.species.q_over_m()

    # V_new = V_old + dt * E
    for i_part=1:self.particle_group.n_particles

       v_new = self%particle_group%get_v(i_part)
       # Evaluate efields at particle position
       xi = self.particle_group.get_x(i_part)
       evaluate(kernel_smoother_1, xi[1], efield_dofs[:,1], efield[1])
       evaluate(kernel_smoother_0, xi[1], efield_dofs[:,2], efield[2])
       v_new = particle_group.get_v(i_part)
       v_new[1:2] = v_new[1:2] + dt * qm * efield
       set_v(particle_group, i_part, v_new)

    end do
    
    # Update bfield
    compute_B_from_E( maxwell_solver, dt, self.efield_dofs[:,2], 
         self.bfield_dofs)
        
end
  
"""
Push H_B: Equations to be solved ``V_new = V_old``

```math
\\begin{eqnarray}
\\partial_t E_1 = 0 & - &> E_{1,new} = E_{1,old} \\
\\partial_t E_2 = - \\partial_{x_1} B & -> & E_{2,new} = E_{2,old}-dt*\\partial_{x_1} B \\
\\partial_t B = 0 & -> & B_new = B_old \\
\\end{eqnarray}
```
"""
function operatorHB(self, dt)
    # Update efield2
    maxwell_solver.compute_E_from_B( dt, 
                                     self.bfield_dofs, 
                                     self.efield_dofs[:,2])
      
end


export HamiltonianSplittingBoris

"""
    HamiltonianSplittingBoris( maxwell_solver,
                               kernel_smoother_0, kernel_smoother_1,
                               particle_group,
                               e_dofs_1, e_dofs_2, b_dofs,
                               x_min, Lx ) 

Boris pusher in GEMPIC framework (spline finite elements)

- `mid` describing the magnetic field at time ``t_{n+1/2}`` (used for push)
- `j_dofs` for kernel representation of current density. 
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
     kernel_smoother_0 :: ParticleMeshCoupling #(order p+1)
     kernel_smoother_1 :: ParticleMeshCoupling #(order p)
     particle_group    :: ParticleGroup
     spline_degree     :: Int64
     Lx                :: Float64
     x_min             :: Float64 # Lower bound for x domain
     delta_x           :: Float64 # Grid spacing

     cell_integrals_0  :: SVector
     cell_integrals_1  :: SVector

     e_dofs     :: Array{Array{Float64, 1},1}
     e_dofs_mid :: Array{Array{Float64, 1},1}
     b_dofs     :: Array{Float64, 1}
     b_dofs_mid :: Array{Float64, 1}
     j_dofs     :: Array{Array{Float64, 1}}

     function HamiltonianSplittingBoris( 
         maxwell_solver :: Maxwell1DFEM,
         kernel_smoother_0 :: ParticleMeshCoupling,
         kernel_smoother_1 :: ParticleMeshCoupling,
         particle_group :: ParticleGroup,
         e_dofs :: Array{Vector{Float64},1},
         b_dofs :: Vector{Float64},
         domain :: Vector{Float64} ) 

         e_dofs_mid = [zeros(Float64, kernel_smoother_1.n_dofs),
                       zeros(Float64, kernel_smoother_1.n_dofs)]

         j_dofs     = [zeros(Float64, kernel_smoother_0.n_dofs),
                       zeros(Float64, kernel_smoother_0.n_dofs)]

         b_dofs_mid = zeros(Float64, kernel_smoother_1.n_dofs)

         x_min         = domain[1]
         Lx            = domain[3]
         spline_degree = 3
         delta_x       = Lx/kernel_smoother_1.n_dofs
    
         cell_integrals_1 = SVector{3}([0.5,  2.0,  0.5] ./ 3.0)
         cell_integrals_0 = SVector{4}([1.0, 11.0, 11.0, 1.0] ./ 24.0)

         new( maxwell_solver, kernel_smoother_0, kernel_smoother_1,
              particle_group, spline_degree, Lx, x_min, delta_x,
              cell_integrals_0, cell_integrals_1,
              e_dofs, e_dofs_mid, b_dofs, b_dofs_mid, j_dofs)

    end

end

export staggering!
"""
    staggering(splitting, dt)

Propagate ``E_0`` to ``E_{1/2}`` and ``x_0`` to ``x_{1/2}`` to initialize 
the staggering
- `splitting` : time splitting object 
- `dt`   : time step
"""
function staggering!(splitting, dt)

    push_x_accumulate_j!(splitting, dt*0.5)

    # (4) Compute E_{n+1}
    splitting.e_dofs_mid[1] .= splitting.e_dofs[1]
    splitting.e_dofs_mid[2] .= splitting.e_dofs[2]
    splitting.j_dofs[1]  .= 0.5dt .* splitting.j_dofs[1]
    splitting.j_dofs[2]  .= 0.5dt .* splitting.j_dofs[2]
    # Ex part:

    compute_e_from_j!(splitting.e_dofs_mid[1], 
                      splitting.maxwell_solver, 
                      splitting.j_dofs[1], 1) 

    #todo "Combine the two steps for efficiency"

    compute_e_from_b!(splitting.e_dofs_mid[2], 
                      splitting.maxwell_solver, 
                      0.5dt, splitting.b_dofs)

    compute_e_from_j!(splitting.e_dofs_mid[2], 
                      splitting.maxwell_solver,
                      splitting.j_dofs[2], 2) 

end


"""
    strang_splitting(splitting, dt, number_steps)

Second order Boris pusher using staggered grid
- `splitting` : time splitting object 
- `dt`   : time step
- `number_steps` : number of time steps
"""
function strang_splitting!(splitting    :: HamiltonianSplittingBoris, 
                           dt           :: Float64, 
                           number_steps :: Int64)

    for i_step = 1:number_steps

        # (1) Compute B_{n+1/2} from B_n
        splitting.b_dofs_mid .= splitting.b_dofs

        compute_b_from_e!(splitting.b_dofs,
                          splitting.maxwell_solver,
                          dt, 
                          splitting.e_dofs_mid[2]) 

        splitting.b_dofs_mid .= (splitting.b_dofs_mid .+ splitting.b_dofs) * .5
       
        # (2) Propagate v: v_{n-1/2} -> v_{n+1/2}
        # (2a) Half time step with E-part
        push_v_epart!( splitting, 0.5dt )
        # (2b) Full time step with B-part
        push_v_bpart!( splitting, dt )
        # (2c) Half time step with E-part
        push_v_epart!( splitting, 0.5dt )
        
        # (3) Propagate x: x_n -> x_{n+1}. Includes also accumulation of j_x, j_y
        push_x_accumulate_j!( splitting, dt )
        
        # (4) Compute E_{n+1}       
        splitting.e_dofs[1] .= splitting.e_dofs_mid[1]
        splitting.e_dofs[2] .= splitting.e_dofs_mid[2]
        splitting.j_dofs[1] .= dt * splitting.j_dofs[1]
        splitting.j_dofs[2] .= dt * splitting.j_dofs[2]

        # Ex part:
        compute_e_from_j!(splitting.e_dofs_mid[1], 
                          splitting.maxwell_solver, 
                          splitting.j_dofs[1], 1)

        # Ey part:    
        compute_e_from_b!(splitting.e_dofs_mid[2], 
                          splitting.maxwell_solver, dt, 
                          splitting.b_dofs) 

        compute_e_from_j!(splitting.e_dofs_mid[2], 
                          splitting.maxwell_solver, 
                          splitting.j_dofs[2], 2)
       
    end
 

end

"""
    push_v_epart(splitting, dt)


Pusher for ``E \\nabla_v `` part

```math
V_{new} = V_{old} + dt ⋅ E
```
"""
function push_v_epart!(splitting, dt)

    qm = splitting.particle_group.q_over_m

    for i_part = 1:splitting.particle_group.n_particles

        xi = get_x(splitting.particle_group, i_part)

        efield_1 = evaluate(splitting.kernel_smoother_1, xi[1], splitting.e_dofs_mid[1])

        efield_2 = evaluate(splitting.kernel_smoother_0, xi[1], splitting.e_dofs_mid[2])

        v_new  = get_v(splitting.particle_group, i_part)
        v_new .= v_new .+ dt * qm .* [efield_1, efield_2]

        set_v(splitting.particle_group, i_part, v_new)

    end
    

end

"""
    push_v_bpart!(splitting, dt)

  Pusher for vxB part
"""
function push_v_bpart!(splitting :: HamiltonianSplittingBoris, 
                       dt        :: Float64)

    qmdt = splitting.particle_group.q_over_m * 0.5 * dt
    v_new = zeros(3)
    
    for i_part=1:splitting.particle_group.n_particles

        vi = get_v(splitting.particle_group, i_part)
        xi = get_x(splitting.particle_group, i_part)

        bfield = evaluate(splitting.kernel_smoother_1, xi[1], 
                          splitting.b_dofs_mid)

        bfield = qmdt * bfield

        M11    = 1.0/(1.0 + bfield^2) 
        M12    = M11 * bfield * 2.0
        M11    = M11 * (1-bfield^2)

        v_new[1] =   M11 * vi[1] + M12 * vi[2]
        v_new[2] = - M12 * vi[1] + M11 * vi[2]
        v_new[3] = 0.0

        set_v(splitting.particle_group, i_part, v_new)

    end

end

"""
    push_x_accumulate_j!(splitting, dt)

Pusher for x and accumulate current densities

For each particle compute the index of the first DoF on the grid it 
contributes to and its position (normalized to cell size one). 

Then update particle position:  ``X_{new} = X_{old} + dt ⋅ V``

!!! note
    `j_dofs` does not hold the values for j itself but 
    for the integrated j.

"""
function push_x_accumulate_j!(splitting, dt)

    n_cells = splitting.kernel_smoother_0.n_dofs

    fill!(splitting.j_dofs[1], 0.0)
    fill!(splitting.j_dofs[2], 0.0)

    for i_part=1:splitting.particle_group.n_particles  

       # Read out particle position and velocity
       x_old = get_x(splitting.particle_group, i_part)
       vi = get_v(splitting.particle_group, i_part)

       # Then update particle position:  X_new = X_old + dt * V
       x_new = x_old .+ dt * vi

       # Get charge for accumulation of j
       wi = get_charge(splitting.particle_group, i_part)
       qoverm = splitting.particle_group.q_over_m

       #todo "Check here also second posibility with sum of two accumulations"
       # Accumulate jx
       add_charge!(splitting.j_dofs[1], splitting.kernel_smoother_1,
                [(x_old[1]+x_new[1])*0.5], wi[1]*vi[1])
       # Accumulate jy
       add_charge!(splitting.j_dofs[2], splitting.kernel_smoother_0,
                [(x_old[1]+x_new[1])*0.5], wi[1]*vi[2])
      
       x_new[1] = mod(x_new[1], splitting.Lx)
       set_x(splitting.particle_group, i_part, x_new)

    end

end 


using StaticArrays

abstract type AbstractSplitting end

export HamiltonianSplitting


"""
    HamiltonianSplitting( maxwell_solver,
                          kernel_smoother_0, kernel_smoother_1,
                          particle_group, e_dofs, b_dofs) 

Hamiltonian splitting type for Vlasov-Maxwell

- Integral over the spline function on each interval (order p+1)
- Integral over the spline function on each interval (order p)
- `e_dofs` describing the two components of the electric field
- `b_dofs` describing the magnetic field
- `j_dofs` for kernel representation of current density. 
"""
struct HamiltonianSplitting

    maxwell_solver    :: AbstractMaxwellSolver
    kernel_smoother_0 :: ParticleMeshCoupling
    kernel_smoother_1 :: ParticleMeshCoupling
    particle_group    :: ParticleGroup

    spline_degree     :: Int
    Lx                :: Float64
    x_min             :: Float64
    delta_x           :: Float64

    cell_integrals_0  :: SVector
    cell_integrals_1  :: SVector

    e_dofs :: Array{Array{Float64,1}}
    b_dofs :: Array{Float64,1}
    j_dofs :: Array{Array{Float64,1}}

    function HamiltonianSplitting( maxwell_solver,
                                   kernel_smoother_0,
                                   kernel_smoother_1,
                                   particle_group,
                                   e_dofs,
                                   b_dofs) 

        # Check that n_dofs is the same for both kernel smoothers.
        @assert kernel_smoother_0.n_dofs == kernel_smoother_1.n_dofs

        j_dofs = [zeros(Float64,kernel_smoother_0.n_dofs) for i in 1:2]

        x_min = maxwell_solver.xmin
        Lx    = maxwell_solver.Lx
        spline_degree = 3
        delta_x = Lx/kernel_smoother_1.n_dofs
        
        cell_integrals_1 = SVector{3}([0.5, 2.0, 0.5] ./ 3.0)
        cell_integrals_0 = SVector{4}([1.0,11.0,11.0,1.0] ./ 24.0)

        new( maxwell_solver, kernel_smoother_0, kernel_smoother_1, 
             particle_group, spline_degree, Lx, x_min, delta_x,
             cell_integrals_0, cell_integrals_1, e_dofs, b_dofs, j_dofs)

    end

end

export strang_splitting!

"""
    strang_splitting( h, dt, number_steps)

Strang splitting
- time splitting object 
- time step
- number of time steps
"""
function strang_splitting!( h            :: HamiltonianSplitting,
                            dt           :: Float64, 
                            number_steps :: Int)

    for i_step = 1:number_steps
       operatorHB(  h, 0.5dt)
       operatorHE(  h, 0.5dt)
       operatorHp2( h, 0.5dt)
       operatorHp1( h, 1.0dt)
       operatorHp2( h, 0.5dt)
       operatorHE(  h, 0.5dt)
       operatorHB(  h, 0.5dt)
    end

end 

export lie_splitting!
"""
    lie_splitting( h, dt, number_steps)

Lie splitting
"""
function lie_splitting!(h            :: HamiltonianSplitting,
                        dt           :: Float64, 
                        number_steps :: Int)

    for i_step = 1:number_steps
       operatorHE(dt)
       operatorHB(dt)
       operatorHp1(dt)
       operatorHp2(dt)
    end

end 

export lie_splitting_back!
"""
    lie_splitting_back(h, dt, number_steps)

Lie splitting (oposite ordering)
"""
function lie_splitting_back!(h            :: HamiltonianSplitting,
                             dt           :: Float64, 
                             number_steps :: Int)

    for i_step = 1:number_steps

       operatorHp2(dt)
       operatorHp1(dt)
       operatorHB(dt)
       operatorHE(dt)

    end

end
  

"""

    operatorHp1(h, dt)

```math
\\begin{aligned}
\\partial_t f + v_1 \\partial_{x_1} f = 0 & \\rightarrow  X_{new} = X_{old} + dt V_1 \\\\
V_{new},2 = V_{old},2 + \\int_0 h V_{old},1 B_{old} & \\\\
\\partial_t E_1 = - \\int v_1 f(t,x_1, v) dv & \\rightarrow E_{1,new} = E_{1,old} - 
\\int \\int v_1 f(t,x_1+s v_1,v) dv ds  & \\\\
\\partial_t E_2 = 0 & \\rightarrow E_{2,new} = E_{2,old} \\\\
\\partial_t B = 0 & \\rightarrow B_{new} = B_{old} 
\\end{aligned}
```

Here we have to accumulate j and integrate over the time interval.
At each ``k=1,...,n_{grid}``, we have for ``s \\in [0,dt]:
j_k(s) =  \\sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h) v_k``,
where ``h`` is the grid spacing and ``N`` the normalized B-spline
 In order to accumulate the integrated ``j``, we normalize the values of 
``x`` to the grid spacing, 
    calling them ``y``, we have
```math
 j_k(s) = \\sum_{i=1,..,N_p} q_i N(y_k+\\frac{s}{h} v_{1,k}-y_i) v_k.
```
 Now, we want the integral

```math
\\int_{0..dt} j_k(s) d s = \\sum_{i=1}^{N_p} q_i v_k 
\\int_{0}{dt} N(y_k+\\frac{s}{h} v_{1,k}-y_i) ds 
=  \\sum_{i=1}^{N_p} q_i v_k  \\int_{0}^{dt}  N(y_k + w v_{1,k}-y_i) dw
```

For each particle compute the index of the first DoF on the grid it contributes 
to and its position (normalized to cell size one). Note: `j_dofs` 
does not hold the values for `j` itself but for the integrated `j`.

Then update particle position:  
``X_{new} = X_{old} + dt * V``
"""
function operatorHp1(h :: HamiltonianSplitting, dt :: Float64)

    n_cells = h.kernel_smoother_0.n_dofs

    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    @inbounds for i_part = 1:h.particle_group.n_particles  

       # Read out particle position and velocity
       x_old = get_x(h.particle_group, i_part)[1]
       vi    = get_v(h.particle_group, i_part)

       # Then update particle position:  X_new = X_old + dt * V
       x_new = x_old + dt * vi[1]

       # Get charge for accumulation of j
       wi     = get_charge(h.particle_group, i_part)
       qoverm = h.particle_group.q_over_m

       vi[2] = add_current_update_v!( h.j_dofs[1], 
                              h.kernel_smoother_1,
                              x_old, 
                              x_new, 
                              wi,
                              qoverm, 
                              h.b_dofs, 
                              vi[2])

       # Accumulate rho for Poisson diagnostics
       add_charge!( h.j_dofs[2],
                    h.kernel_smoother_0, 
                    x_new, 
                    wi)
      
       x_new = mod(x_new, h.Lx)

       set_x(h.particle_group, i_part, x_new)
       set_v(h.particle_group, i_part, vi)

    end

    # Update the electric field.
    compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.j_dofs[1], 1)

end

"""
    operatorHp2(h, dt)

Push Hp2: Equations to solve are

```math
\\begin{aligned} X_{new}  =  X_{old} & \\\\
V_{new,1} = V_{old,1} + \\int_0 h V_{old,2} B_{old} & \\\\
\\partial_t E_1 = 0 & \\rightarrow  E_{1,new} = E_{1,old}  \\\\
\\partial_t E_2 = - \\int v_2 f(t,x_1, v) dv & \\rightarrow 
E_{2,new} = E_{2,old} - \\int \\int v_2 f(t,x_1 + s v_1,v) dv ds\\\\
\\partial_t B = 0 & => & B_{new} = B_{old} & \\\\
\\end{aligned}
```
"""
function operatorHp2(h :: HamiltonianSplitting, dt :: Float64)
    
    n_cells = h.kernel_smoother_0.n_dofs

    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    qm = h.particle_group.q_over_m

    # Update v_1
    for i_part=1:h.particle_group.n_particles

       # Evaluate b at particle position (splines of order p)
       xi    = get_x(h.particle_group, i_part)[1]
       b     = evaluate(h.kernel_smoother_1, xi[1], h.b_dofs)
       vi    = get_v(h.particle_group, i_part)
       vi[1] = vi[1] + dt * qm * vi[2] * b
       set_v(h.particle_group, i_part, vi)

       xi = get_x(h.particle_group, i_part)[1]

       # Scale vi by weight to combine both factors 
       #for accumulation of integral over j

       wi = get_charge(h.particle_group, i_part) * vi[2]

       add_charge!( h.j_dofs[2], h.kernel_smoother_0, xi, wi[1])

    end

    # Update the electric field. Also, we still need to scale with 1/Lx 
    
    h.j_dofs[2] .= h.j_dofs[2] .* dt

    compute_e_from_j!( h.e_dofs[2], h.maxwell_solver, h.j_dofs[2], 2)
    
end
  
"""
    operatorHE(h, dt)

Push ``H_E``: Equations to be solved
```math
\\begin{aligned}
\\partial_t f + E_1 \\partial_{v_1} f + E_2 \\partial_{v_2} f = 0 &\\rightarrow& V_{new} = V_{old} + dt * E \\\\
\\partial_t E_1 = 0 &\\rightarrow& E_{1,new} = E_{1,old} \\\\
\\partial_t E_2 = 0 &\\rightarrow& E_{2,new} = E_{2,old} \\\\
\\partial_t B + \\partial_{x_1} E_2 = 0 &\\rightarrow& B_{new} = B_{old} - dt \\partial_{x_1} E_2
\\end{aligned}
```
"""
function operatorHE(h :: HamiltonianSplitting, dt :: Float64)

    qm = h.particle_group.q_over_m

    # V_new = V_old + dt * E
    for i_part=1:h.particle_group.n_particles

       v_new = get_v( h.particle_group, i_part)
       # Evaluate efields at particle position
       xi = get_x(h.particle_group, i_part)
       e1 = evaluate(h.kernel_smoother_1, xi[1], h.e_dofs[1])
       e2 = evaluate(h.kernel_smoother_0, xi[1], h.e_dofs[2])
       v_new = get_v(h.particle_group, i_part)
       v_new[1:2] .= v_new[1:2] .+ dt * qm * [e1, e2]
       set_v(h.particle_group, i_part, v_new)

    end
    
    # Update bfield
    compute_b_from_e!( h.b_dofs, h.maxwell_solver, dt, h.e_dofs[2])
        
end
  
@doc raw"""
    operatorHB(h, dt)

Push ``H_B``: Equations to be solved ``V_{new} = V_{old}``

```math
\begin{aligned}
\partial_t E_1 = 0 & \rightarrow & E_{1,new} = E_{1,old} \\
\partial_t E_2 = - \partial_{x_1} B & \rightarrow & E_{2,new} = E_{2,old}-dt*\partial_{x_1} B \\
\partial_t B = 0 & \rightarrow & B_{new} = B_{old} \\
\end{aligned}
```
"""
function operatorHB(h :: HamiltonianSplitting, dt :: Float64)
    compute_e_from_b!( h.e_dofs[2], h.maxwell_solver, dt, h.b_dofs)
end

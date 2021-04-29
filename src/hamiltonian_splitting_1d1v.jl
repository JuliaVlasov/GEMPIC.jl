export strang_splitting!

"""
    strang_splitting( h, dt, number_steps)

Strang splitting
- time splitting object 
- time step
- number of time steps
"""
function strang_splitting!( h            :: HamiltonianSplitting{1,1},
                            dt           :: Float64, 
                            number_steps :: Int64)

    for i_step = 1:number_steps
       operatorHB(  h, 0.5dt)
       operatorHp1( h, dt)
       operatorHB(  h, 0.5dt)
    end

end 

"""

    operatorHp1(h, dt)

```math
\\begin{eqnarray}
\\partial_t f + v_1 \\partial_{x_1} f = 0 & -> &  X_{new} = X_{old} + dt V_1 \\\\
V_{new},2 = V_{old},2 + \\int_0 h V_{old},1 B_{old} && \\\\
\\partial_t E_1 = - \\int v_1 f(t,x_1, v) dv & -> & E_{1,new} = E_{1,old} - 
\\int \\int v_1 f(t,x_1+s v_1,v) dv ds  && \\\\
\\partial_t E_2 = 0 & -> & E_{2,new} = E_{2,old} \\\\
\\partial_t B = 0 & => & B_{new} = B_{old} 
\\end{eqnarray}
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
function operatorHp1(h :: HamiltonianSplitting{1,1}, dt :: Float64)

    n_cells = h.kernel_smoother_0.n_dofs

    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    for i_part = 1:h.particle_group.n_particles  

       # Read out particle position and velocity
       x_old = h.particle_group.array[1, i_part]
       v_old = h.particle_group.array[2, i_part]
        
       # Then update particle position:  X_new = X_old + dt * V
       x_new = x_old + dt * v_old

       # Get charge for accumulation of j
       wi     = get_charge(h.particle_group, i_part)
       qoverm = h.particle_group.q_over_m

       add_current_update_v!( h.j_dofs[1], 
                              h.kernel_smoother_1,
                              x_old, 
                              x_new, 
                              wi[1],
                              qoverm, 
                              v_old)

       x_new = mod(x_new, h.Lx)
       h.particle_group.array[1, i_part] = x_new
       
    end

    compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.j_dofs[1], 1)

end

"""
    operatorHB(h, dt)

Push H_B: Equations to be solved ``V_{new} = V_{old}``

```math
\\begin{eqnarray}
\\partial_t E_1 = 0 & -> & E_{1,new} = E_{1,old} \\\\
\\partial_t E_2 = - \\partial_{x_1} B & -> & E_{2,new} = E_{2,old}-dt*\\partial_{x_1} B \\\\
\\partial_t B = 0 & -> & B_{new} = B_{old} \\\\
\\end{eqnarray}
```
"""
function operatorHB(h :: HamiltonianSplitting{1,1}, dt :: Float64)

    for i_part=1:h.particle_group.n_particles

        # Evaluate b at particle position (splines of order p)
        xi = h.particle_group.array[1, i_part]
        vi = h.particle_group.array[2, i_part]
        e1 = evaluate(h.kernel_smoother_1, xi, h.e_dofs[1])
        vi = vi + dt  * e1;  
        h.particle_group.array[2,i_part] = vi
    end

end

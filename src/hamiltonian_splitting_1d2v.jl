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
function operatorHp1(h::HamiltonianSplitting{1,2}, dt::Float64)
    np = h.particle_group.n_particles
    n_cells = h.kernel_smoother_0.n_dofs

    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    tasks = map(h.chunks) do i_chunk
        @spawn begin
            buffer = zero(h.j_dofs[1])
            for i_part in i_chunk

                # Read out particle position and velocity
                x_old = h.particle_group.array[1, i_part]
                v1_old = h.particle_group.array[2, i_part]
                v2_old = h.particle_group.array[3, i_part]

                # Then update particle position:  X_new = X_old + dt * V
                x_new = x_old + dt * v1_old

                # Get charge for accumulation of j
                wi = h.particle_group.array[4, i_part]
                wi = wi * h.particle_group.charge
                wi = wi * h.particle_group.common_weight

                qoverm = h.particle_group.q_over_m

                v2_new = add_current_update_v!(
                    buffer,
                    h.kernel_smoother_1,
                    x_old,
                    x_new,
                    wi,
                    qoverm,
                    h.b_dofs,
                    v2_old,
                )

                x_new = mod(x_new, h.Lx)

                h.particle_group.array[1, i_part] = x_new
                h.particle_group.array[3, i_part] = v2_new
            end
            return buffer
        end
    end

    h.j_dofs[1] .= reduce(+, fetch.(tasks))

    tasks = map(h.chunks) do i_chunk 
        @spawn begin
            buffer = zero(h.j_dofs[2])
            for i_part in i_chunk

                # Read out particle position and charge
                x = h.particle_group.array[1, i_part]
                w = h.particle_group.array[4, i_part]
                w = w * h.particle_group.charge
                w = w * h.particle_group.common_weight

                # Accumulate rho for Poisson diagnostics
                add_charge!(buffer, h.kernel_smoother_0, x, w)
            end
            return buffer
        end
    end

    h.j_dofs[2] .= reduce(+, fetch.(tasks))

    # Update the electric field.
    return compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.j_dofs[1], 1)
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
function operatorHp2(h::HamiltonianSplitting{1,2}, dt::Float64)
    n_cells = h.kernel_smoother_0.n_dofs

    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    qm = h.particle_group.q_over_m

    np = h.particle_group.n_particles

    # Update v_1

    tasks = map(h.chunks) do i_chunk 
        @spawn begin
            buffer = zero(h.j_dofs[2])
            for i_part in i_chunk

                # Evaluate b at particle position (splines of order p)
                x1 = h.particle_group.array[1, i_part]
                v1 = h.particle_group.array[2, i_part]
                v2 = h.particle_group.array[3, i_part]

                b = evaluate(h.kernel_smoother_1, x1, h.b_dofs)
                v1 = v1 + dt * qm * v2 * b

                h.particle_group.array[2, i_part] = v1

                # Scale vi by weight to combine both factors 
                #for accumulation of integral over j
                w = h.particle_group.array[4, i_part]
                w = w * h.particle_group.charge
                w = w * h.particle_group.common_weight
                w = w * v2

                add_charge!(buffer, h.kernel_smoother_0, x1, w)
            end
            return buffer
        end
    end

    h.j_dofs[2] .= reduce(+, fetch.(tasks))

    # Update the electric field. Also, we still need to scale with 1/Lx 

    h.j_dofs[2] .= h.j_dofs[2] .* dt

    return compute_e_from_j!(h.e_dofs[2], h.maxwell_solver, h.j_dofs[2], 2)
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
function operatorHE(h::HamiltonianSplitting{1,2}, dt::Float64)
    qm = h.particle_group.q_over_m

    np = h.particle_group.n_particles

    # V_new = V_old + dt * E

    @threads for i_part in 1:np

        v_old1 = h.particle_group.array[2, i_part]
        v_old2 = h.particle_group.array[3, i_part]

        # Evaluate efields at particle position

        xi = h.particle_group.array[1, i_part]

        e1 = evaluate(h.kernel_smoother_1, xi, h.e_dofs[1])
        e2 = evaluate(h.kernel_smoother_0, xi, h.e_dofs[2])

        v_new1 = v_old1 + dt * qm * e1
        v_new2 = v_old2 + dt * qm * e2

        h.particle_group.array[2, i_part] = v_new1
        h.particle_group.array[3, i_part] = v_new2
    end

    # Update bfield
    return compute_b_from_e!(h.b_dofs, h.maxwell_solver, dt, h.e_dofs[2])
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
function operatorHB(h::HamiltonianSplitting, dt::Float64)
    return compute_e_from_b!(h.e_dofs[2], h.maxwell_solver, dt, h.b_dofs)
end

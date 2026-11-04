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
    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    if length(h.hp1_workspaces) == 1
        deposit_hp1!(h, dt, first(h.chunks), h.hp1_workspaces[1])
    else
        @sync for (i, i_chunk) in enumerate(h.chunks)
            @spawn deposit_hp1!(h, dt, i_chunk, h.hp1_workspaces[i])
        end
    end

    for workspace in h.hp1_workspaces
        h.j_dofs[1] .+= workspace.current
        h.j_dofs[2] .+= workspace.charge
    end

    return compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.j_dofs[1], 1)
end

function deposit_hp1!(h::HamiltonianSplitting{1,2}, dt, i_chunk, workspace::Hp1Workspace)
    fill!(workspace.current, 0.0)
    fill!(workspace.charge, 0.0)
    particles = h.particle_group.array
    @inbounds for i_part in i_chunk
        x_old = particles[1, i_part]
        v1_old = particles[2, i_part]
        v2_old = particles[3, i_part]
        x_new = x_old + dt * v1_old

        wi = particles[4, i_part] * h.particle_group.charge
        wi *= h.particle_group.common_weight

        v2_new = add_current_update_v!(
            workspace.current,
            workspace.kernel_1,
            x_old,
            x_new,
            wi,
            h.particle_group.q_over_m,
            h.b_dofs,
            v2_old,
        )

        x_new = mod(x_new, h.Lx)
        particles[1, i_part] = x_new
        particles[3, i_part] = v2_new

        # Accumulate rho for Poisson diagnostics at the updated position.
        add_charge!(workspace.charge, workspace.kernel_0, x_new, wi)
    end
    return nothing
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
    # Convert once per field, then share read-only coefficients between particle tasks.
    b_to_pp!(h.he_pp[1], h.kernel_smoother_1.spline_pp, h.kernel_smoother_1.n_grid, h.e_dofs[1])
    b_to_pp!(h.he_pp[2], h.kernel_smoother_0.spline_pp, h.kernel_smoother_0.n_grid, h.e_dofs[2])

    np = h.particle_group.n_particles
    n_chunks = length(h.hp1_workspaces)
    if n_chunks == 1
        push_he!(h, dt, 1:np)
    else
        chunk_size = np ÷ n_chunks
        @threads for i in 1:n_chunks
            push_he!(h, dt, ((i - 1) * chunk_size + 1):(i * chunk_size))
        end
    end

    # Update bfield
    return compute_b_from_e!(h.b_dofs, h.maxwell_solver, dt, h.e_dofs[2])
end

function push_he!(h::HamiltonianSplitting{1,2}, dt, i_chunk)
    particles = h.particle_group.array
    qm = h.particle_group.q_over_m
    e1_pp, e2_pp = h.he_pp
    @inbounds for i_part in i_chunk
        xi = particles[1, i_part]
        e1 = evaluate_pp(h.kernel_smoother_1, xi, e1_pp)
        e2 = evaluate_pp(h.kernel_smoother_0, xi, e2_pp)
        particles[2, i_part] += dt * qm * e1
        particles[3, i_part] += dt * qm * e2
    end
    return nothing
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

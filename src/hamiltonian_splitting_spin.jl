export HamiltonianSplittingSpin

"""
    HamiltonianSplittingSpin( maxwell_solver,
                              kernel_smoother_0, kernel_smoother_1,
                              particle_group, e_dofs, a_dofs ) 

Hamiltonian splitting type for Vlasov-Maxwell

- Integral over the spline function on each interval (order p+1)
- Integral over the spline function on each interval (order p)
- `e_dofs` describing the two components of the electric field
- `a_dofs` describing the potential vector
"""
struct HamiltonianSplittingSpin
    maxwell_solver::AbstractMaxwellSolver
    kernel_smoother_0::ParticleMeshCoupling1D
    kernel_smoother_1::ParticleMeshCoupling1D
    kernel_smoother_2::ParticleMeshCoupling1D
    particle_group::ParticleGroup

    spline_degree::Int64
    Lx::Float64
    x_min::Float64
    delta_x::Float64

    cell_integrals_0::SVector
    cell_integrals_1::SVector

    e_dofs::Array{Array{Float64,1}}
    a_dofs::Array{Array{Float64,1}}
    j_dofs::Array{Array{Float64,1}}
    part::Array{Array{Float64,1}}

    HH::Float64

    function HamiltonianSplittingSpin(
        maxwell_solver,
        kernel_smoother_0,
        kernel_smoother_1,
        kernel_smoother_2,
        particle_group,
        e_dofs,
        a_dofs,
        HH=0.0,
    ) #0.00022980575) 

        # Check that n_dofs is the same for both kernel smoothers.
        @assert kernel_smoother_0.n_dofs == kernel_smoother_1.n_dofs

        j_dofs = [zeros(Float64, kernel_smoother_0.n_dofs) for i in 1:2]

        nx = maxwell_solver.n_dofs
        part = [zeros(Float64, nx) for i in 1:4]
        x_min = maxwell_solver.xmin
        Lx = maxwell_solver.Lx
        spline_degree = 3
        delta_x = Lx / kernel_smoother_1.n_dofs

        cell_integrals_1 = SVector{3}([0.5, 2.0, 0.5] ./ 3.0)
        cell_integrals_0 = SVector{4}([1.0, 11.0, 11.0, 1.0] ./ 24.0)

        return new(
            maxwell_solver,
            kernel_smoother_0,
            kernel_smoother_1,
            kernel_smoother_2,
            particle_group,
            spline_degree,
            Lx,
            x_min,
            delta_x,
            cell_integrals_0,
            cell_integrals_1,
            e_dofs,
            a_dofs,
            j_dofs,
            part,
            HH,
        )
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
function strang_splitting!(h::HamiltonianSplittingSpin, dt::Float64, number_steps::Int64)
    for i_step in 1:number_steps
        operatorHE(h, 0.5dt)
        operatorHp(h, 0.5dt)
        operatorHA(h, 0.5dt)
        operatorHs(h, dt)
        operatorHA(h, 0.5dt)
        operatorHp(h, 0.5dt)
        operatorHE(h, 0.5dt)
    end
end

"""

    operatorHp(h, dt)

```math
\\begin{aligned}
\\dot{x} & =p \\\\
\\dot{E}_x & = - \\int (p f ) dp ds
\\end{aligned}
```
"""
function operatorHp(h::HamiltonianSplittingSpin, dt::Float64)
    n_cells = h.kernel_smoother_0.n_dofs

    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    for i_part in 1:(h.particle_group.n_particles)

        # Read out particle position and velocity
        x_old = h.particle_group.array[1, i_part]
        v_old = h.particle_group.array[2, i_part]

        # Then update particle position:  X_new = X_old + dt * V
        x_new = x_old + dt * v_old

        # Get charge for accumulation of j
        w = get_charge(h.particle_group, i_part)

        add_current_update_v!(h.j_dofs[1], h.kernel_smoother_1, x_old, x_new, w)

        x_new = mod(x_new, h.Lx)
        h.particle_group.array[1, i_part] = x_new
    end

    # Update the electric field.
    return compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.j_dofs[1], 1)
end

"""
    operatorHA(h, dt)

```math
\\begin{aligned}
\\dot{p} = (A_y, A_z) \\cdot \\partial_x (A_y, A_z)   \\\\
\\dot{Ey} = -\\partial_x^2 A_y + A_y \\rho \\\\
\\dot{Ez} = -\\partial_x^2 A_z + A_z \\rho \\\\
\\end{aligned}
```
"""
function operatorHA(h::HamiltonianSplittingSpin, dt::Float64)
    n_cells = h.kernel_smoother_0.n_dofs

    fill!.(h.part, 0.0)

    qm = h.particle_group.q_over_m
    # Update velocities
    aa = zeros(Float64, n_cells)

    tmp_jdofs1 = zeros(h.kernel_smoother_1.n_span + 1)
    tmp_jdofs2 = zeros(h.kernel_smoother_0.n_span)

    for i_part in 1:(h.particle_group.n_particles)

        # fill!(h.j_dofs[1], 0.0)
        # fill!(h.j_dofs[2], 0.0)

        # Evaluate b at particle position (splines of order p)
        xi = h.particle_group.array[1, i_part]
        vi = h.particle_group.array[2, i_part]
        wi = get_charge(h.particle_group, i_part)

        # add_charge!( h.j_dofs[2], h.kernel_smoother_0, xi, 1.0)

        xn = (xi - h.kernel_smoother_0.xmin) / h.kernel_smoother_0.delta_x
        index = trunc(Int, xn)
        xn = xn - index
        index = index - h.kernel_smoother_0.spline_degree

        uniform_bsplines_eval_basis!(
            h.kernel_smoother_0.spline_val, h.kernel_smoother_0.spline_degree, xn
        )

        for i in 1:(h.kernel_smoother_0.n_span)
            tmp_jdofs2[i] = h.kernel_smoother_0.spline_val[i] * h.kernel_smoother_0.scaling
        end

        # add_charge!( h.j_dofs[1], h.kernel_smoother_1, xi, 1.0)

        xn = (xi - h.kernel_smoother_1.xmin) / h.kernel_smoother_1.delta_x
        index = trunc(Int, xn)
        xn = xn - index
        index = index - h.kernel_smoother_1.spline_degree

        uniform_bsplines_eval_basis!(
            h.kernel_smoother_1.spline_val, h.kernel_smoother_1.spline_degree, xn
        )

        for i in 1:(h.kernel_smoother_1.n_span)
            tmp_jdofs1[i] = h.kernel_smoother_1.spline_val[i] * h.kernel_smoother_1.scaling
        end

        # values of the derivatives of basis function
        # compute_rderivatives_from_basis!(aa, h.maxwell_solver, h.j_dofs[1])

        coef = 1 / h.maxwell_solver.delta_x
        # relation betwen spline coefficients for strong Ampere
        for i in 1:(h.kernel_smoother_1.n_span)
            index1d = mod1(index + i, n_cells)
            aa[index1d] = coef * (tmp_jdofs1[i] - tmp_jdofs1[i + 1])
        end
        ## treat Periodic point
        #aa[end] =  coef * ( field_in[end] - field_in[1] )
        # h.j_dofs[1] .= aa

        for i in 1:(h.kernel_smoother_1.n_span)
            index1d = mod1(index + i, n_cells)
            tmp_jdofs1[i] = aa[index1d]
        end

        #p11 = h.a_dofs[1]'h.j_dofs[1]
        #p21 = h.a_dofs[2]'h.j_dofs[1]

        #p12 = h.a_dofs[1]'h.j_dofs[2]
        #p22 = h.a_dofs[2]'h.j_dofs[2]

        p11 = 0.0
        p21 = 0.0
        p12 = 0.0
        p22 = 0.0
        for i in 1:(h.kernel_smoother_0.n_span)
            index1d = mod1(index + i, n_cells)
            p11 += h.a_dofs[1][index1d] * tmp_jdofs1[i]
            p21 += h.a_dofs[2][index1d] * tmp_jdofs1[i]
            p12 += h.a_dofs[1][index1d] * tmp_jdofs2[i]
            p22 += h.a_dofs[2][index1d] * tmp_jdofs2[i]
        end

        vi = vi - dt * (p11 * p12 + p22 * p21)

        h.particle_group.array[2, i_part] = vi

        # below we solve electric field
        # first define part1 and part2 to be 0 vector
        for i in 1:(h.kernel_smoother_0.n_span)
            index1d = mod1(index + i, n_cells)
            h.part[1][index1d] = h.part[1][index1d] - tmp_jdofs2[i] * dt * wi * p12
            h.part[2][index1d] = h.part[2][index1d] - tmp_jdofs2[i] * dt * wi * p22
        end
    end

    # Update the electric field. Also, we still need to scale with 1/Lx 

    compute_e_from_j!(h.e_dofs[2], h.maxwell_solver, h.part[1], 2)
    compute_e_from_j!(h.e_dofs[3], h.maxwell_solver, h.part[2], 2)

    # define part3 and part4
    compute_lderivatives_from_basis!(h.part[3], h.maxwell_solver, h.a_dofs[1])
    compute_lderivatives_from_basis!(h.part[4], h.maxwell_solver, h.a_dofs[2])

    compute_e_from_b!(h.e_dofs[2], h.maxwell_solver, dt, h.part[3])
    return compute_e_from_b!(h.e_dofs[3], h.maxwell_solver, dt, h.part[4])
end

"""
    operatorHE(h, dt)

```math
\\begin{aligned}
\\dot{v}   & =  E_x \\\\
\\dot{A}_y & = -E_y \\\\
\\dot{A}_z & = -E_z
\\end{aligned}
```
"""
function operatorHE(h::HamiltonianSplittingSpin, dt::Float64)
    np = h.particle_group.n_particles

    if np > 2
        chunks = Iterators.partition(1:np, np ÷ nthreads())
    else
        chunks = Iterators.partition(1:np, nthreads())
    end

    @sync for i_chunk in chunks
        @spawn begin
            for i_part in i_chunk
                xi = h.particle_group.array[1, i_part]
                vi = h.particle_group.array[2, i_part]

                e1 = evaluate(h.kernel_smoother_1, xi, h.e_dofs[1])

                h.particle_group.array[2, i_part] = vi + dt * e1
            end
        end
    end

    h.a_dofs[1] .-= dt * h.e_dofs[2]
    return h.a_dofs[2] .-= dt * h.e_dofs[3]
end

"""
    operatorHs(h, dt)

Push H_s: Equations to be solved
```math
\\begin{aligned}
\\dot{s} &= s x B = (s_y \\partial_x A_y +s_z \\partial_x Az, -s_x \\partial_x A_y, -s_x \\partial_x A_z)  \\\\
\\dot{p} &= s \\cdot \\partial_x B = -s_y \\partial^2_{x} A_z + s_z \\partial^2_{x} A_y \\\\
\\dot{E}_y &=   \\int (s_z \\partial_x f) dp ds \\\\
\\dot{E}_z &= - \\int (s_y \\partial_x f) dp ds 
\\end{aligned}
```
"""
function operatorHs(h::HamiltonianSplittingSpin, dt::Float64)
    n_cells = h.kernel_smoother_0.n_dofs

    fill!.(h.part, 0.0)

    hat_v = zeros(Float64, 3, 3)

    S = zeros(Float64, 3)
    St = zeros(Float64, 3)
    V = zeros(Float64, 3)
    aa = zeros(Float64, n_cells)

    @inbounds for i_part in 1:(h.particle_group.n_particles)
        xi = h.particle_group.array[1, i_part]
        vi = h.particle_group.array[2, i_part]

        # Evaluate efields at particle position
        fill!(h.j_dofs[1], 0.0)
        fill!(h.j_dofs[2], 0.0)

        add_charge!(h.j_dofs[2], h.kernel_smoother_1, xi, 1.0)

        compute_rderivatives_from_basis!(h.j_dofs[1], h.maxwell_solver, h.j_dofs[2])

        Y = h.a_dofs[1]'h.j_dofs[1]
        Z = h.a_dofs[2]'h.j_dofs[1]
        V .= [0, Z, -Y]

        hat_v[1, 2] = Y
        hat_v[1, 3] = Z
        hat_v[2, 1] = -Y
        hat_v[3, 1] = -Z

        s1 = h.particle_group.array[4, i_part]
        s2 = h.particle_group.array[5, i_part]
        s3 = h.particle_group.array[6, i_part]

        vnorm = norm(V)

        if vnorm > 1e-14
            S .= (
                [s1, s2, s3] .+
                (
                    sin(dt * vnorm) / vnorm * hat_v +
                    0.5 * (sin(dt / 2 * vnorm) / (vnorm / 2))^2 * hat_v^2
                ) * [s1, s2, s3]
            )
        else
            S .= [s1, s2, s3]
        end

        h.particle_group.array[4, i_part] = S[1]
        h.particle_group.array[5, i_part] = S[2]
        h.particle_group.array[6, i_part] = S[3]

        if vnorm > 1e-14
            St .= (
                dt .* [s1, s2, s3] +
                (
                    2 * (sin(dt * vnorm / 2) / vnorm)^2 * hat_v +
                    2.0 / (vnorm^2) * (dt / 2 - sin(dt * vnorm) / 2 / vnorm) * hat_v^2
                ) * [s1, s2, s3]
            )
        else
            St .= dt .* [s1, s2, s3]
        end

        wi = get_charge(h.particle_group, i_part)

        # define part1 and part2
        h.part[1] .+= wi[1] * St[3] * h.j_dofs[2]
        h.part[2] .+= wi[1] * St[2] * h.j_dofs[2]

        # update velocity
        fill!(h.j_dofs[2], 0.0)
        add_charge!(h.j_dofs[2], h.kernel_smoother_2, xi, 1.0)
        compute_rderivatives_from_basis!(h.j_dofs[1], h.maxwell_solver, h.j_dofs[2])
        compute_rderivatives_from_basis!(aa, h.maxwell_solver, h.j_dofs[1])
        h.j_dofs[1] .= aa
        vi =
            vi -
            h.HH * (h.a_dofs[2]' * h.j_dofs[1] * St[2] + h.a_dofs[1]' * h.j_dofs[1] * St[3])

        set_v!(h.particle_group, i_part, vi)
    end

    # Update bfield
    compute_rderivatives_from_basis!(h.j_dofs[1], h.maxwell_solver, h.part[1])
    compute_rderivatives_from_basis!(h.j_dofs[2], h.maxwell_solver, -h.part[2])

    h.j_dofs[1] .*= h.HH
    h.j_dofs[2] .*= h.HH

    compute_e_from_j!(h.e_dofs[2], h.maxwell_solver, h.j_dofs[1], 2)
    return compute_e_from_j!(h.e_dofs[3], h.maxwell_solver, h.j_dofs[2], 2)
end

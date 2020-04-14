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

    maxwell_solver    :: AbstractMaxwellSolver
    kernel_smoother_0 :: ParticleMeshCoupling
    kernel_smoother_1 :: ParticleMeshCoupling
    kernel_smoother_2 :: ParticleMeshCoupling
    particle_group    :: ParticleGroup

    spline_degree     :: Int64
    Lx                :: Float64
    x_min             :: Float64
    delta_x           :: Float64

    cell_integrals_0  :: SVector
    cell_integrals_1  :: SVector

    e_dofs :: Array{Array{Float64,1}}
    a_dofs :: Array{Array{Float64,1}}
    j_dofs :: Array{Array{Float64,1}}
    part   :: Array{Array{Float64,1}}

    function HamiltonianSplittingSpin( maxwell_solver,
                                   kernel_smoother_0,
                                   kernel_smoother_1,
                                   kernel_smoother_2,
                                   particle_group,
                                   e_dofs,
                                   a_dofs) 

        # Check that n_dofs is the same for both kernel smoothers.
        @assert kernel_smoother_0.n_dofs == kernel_smoother_1.n_dofs

        j_dofs = [zeros(Float64,kernel_smoother_0.n_dofs) for i in 1:2]

        nx = maxwell_solver.n_dofs
        part = [zeros(Float64, nx) for i in 1:4]
        x_min = maxwell_solver.xmin
        Lx    = maxwell_solver.Lx
        spline_degree = 3
        delta_x = Lx/kernel_smoother_1.n_dofs
        
        cell_integrals_1 = SVector{3}([0.5, 2.0, 0.5] ./ 3.0)
        cell_integrals_0 = SVector{4}([1.0,11.0,11.0,1.0] ./ 24.0)

        new( maxwell_solver, kernel_smoother_0, kernel_smoother_1, kernel_smoother_2, 
             particle_group, spline_degree, Lx, x_min, delta_x,
             cell_integrals_0, cell_integrals_1, e_dofs, a_dofs, j_dofs,
             part)

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
function strang_splitting!( h           :: HamiltonianSplittingSpin,
                           dt           :: Float64, 
                           number_steps :: Int64)

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

"""

    operatorHp1(h, dt)

```math
\\begin{aligned}
\\partial_t f + v_1 \\partial_{x_1} f = 0 & -> &  X_{new} = X_{old} + dt V_1 \\\\
V_{new},2 = V_{old},2 + \\int_0 h V_{old},1 B_{old} && \\\\
\\partial_t E_1 = - \\int v_1 f(t,x_1, v) dv & -> & E_{1,new} = E_{1,old} - 
\\int \\int v_1 f(t,x_1+s v_1,v) dv ds  && \\\\
\\partial_t E_2 = 0 & -> & E_{2,new} = E_{2,old} \\\\
\\partial_t B = 0 & => & B_{new} = B_{old} 
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
function operatorHp1(h :: HamiltonianSplittingSpin, dt :: Float64)

    n_cells = h.kernel_smoother_0.n_dofs

    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    for i_part = 1:h.particle_group.n_particles  

       # Read out particle position and velocity
       x_old = h.particle_group.particle_array[1, i_part]
       v_old = h.particle_group.particle_array[2, i_part]
        
       # Then update particle position:  X_new = X_old + dt * V
       x_new = x_old + dt * v_old

       # Get charge for accumulation of j
       w = get_charge(h.particle_group, i_part)

       add_current_update_v!( h.j_dofs[1], 
                              h.kernel_smoother_1,
                              x_old, 
                              x_new, 
                              w)

       x_new = mod(x_new, h.Lx)
       h.particle_group.particle_array[1, i_part] = x_new
       
    end

    # Update the electric field.
    compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.j_dofs[1], 1)

end

function operatorHp2(h :: HamiltonianSplittingSpin, dt :: Float64)
    
    n_cells = h.kernel_smoother_0.n_dofs

    fill!.(h.part, 0.0)
    
    qm = h.particle_group.q_over_m
    # Update v_1
    aa = zeros(Float64,n_cells)
    for i_part=1:h.particle_group.n_particles

        fill!(h.j_dofs[1], 0.0)
        fill!(h.j_dofs[2], 0.0)
        # Evaluate b at particle position (splines of order p)
        xi = h.particle_group.particle_array[1, i_part]
        vi = h.particle_group.particle_array[2, i_part]
        wi = get_charge(h.particle_group, i_part) 

        add_charge!( h.j_dofs[2], h.kernel_smoother_0, xi, 1.0)# R0 
        add_charge!( h.j_dofs[1], h.kernel_smoother_1, xi, 1.0)# R1

        # values of the derivatives of basis function
        compute_rderivatives_from_basis!(aa, h.maxwell_solver, h.j_dofs[1])

        h.j_dofs[1] .= aa

        p11 = h.a_dofs[1]'h.j_dofs[1]
        p12 = h.a_dofs[1]'h.j_dofs[2]
        p21 = h.a_dofs[2]'h.j_dofs[1]
        p22 = h.a_dofs[2]'h.j_dofs[2]

        vi = vi - dt * ( p11 * p12 + p22 * p21 )
        
        h.particle_group.particle_array[2, i_part] = vi
        
        # below we solve electric field
        # first define part1 and part2 to be 0 vector

        @. h.j_dofs[1]  = h.j_dofs[2]
        @. h.j_dofs[1] *= dt * wi * p12 
        @. h.part[1]   -= h.j_dofs[1]

        @. h.j_dofs[2] *= dt * wi * p22 
        @. h.part[2]   -= h.j_dofs[2]
        
    end

    # Update the electric field. Also, we still need to scale with 1/Lx 

    compute_e_from_j!( h.e_dofs[2], h.maxwell_solver, h.part[1], 2)
    compute_e_from_j!( h.e_dofs[3], h.maxwell_solver, h.part[2], 2)

    # define part3 and part4
    compute_lderivatives_from_basis!(h.part[3], h.maxwell_solver, h.a_dofs[1])
    compute_lderivatives_from_basis!(h.part[4], h.maxwell_solver, h.a_dofs[2])
    
    compute_e_from_b!( h.e_dofs[2], h.maxwell_solver, dt, h.part[3])
    compute_e_from_b!( h.e_dofs[3], h.maxwell_solver, dt, h.part[4])
    
    
end


"""
    operatorHB(h, dt)

Push ``H_B``: Equations to be solved ``V_{new} = V_{old}``

```math
\\begin{aligned}
\\partial_t E_1 = 0 & -> & E_{1,new} = E_{1,old} \\\\
\\partial_t E_2 = - \\partial_{x_1} B & -> & E_{2,new} = E_{2,old}-dt*\\partial_{x_1} B \\\\
\\partial_t B = 0 & -> & B_{new} = B_{old} \\\\
\\end{aligned}
```
"""
function operatorHB(h :: HamiltonianSplittingSpin, dt :: Float64)

    np = h.particle_group.n_particles

    @sync for i_chunk in Iterators.partition(1:np, nthreads())
        @spawn begin
            for i_part in i_chunk

                xi = h.particle_group.particle_array[1, i_part]
                vi = h.particle_group.particle_array[2, i_part]

                e1 = evaluate(h.kernel_smoother_1, xi, h.e_dofs[1])

                h.particle_group.particle_array[2, i_part] = vi + dt * e1
            end
        end
    end

    h.a_dofs[1] .-= dt*h.e_dofs[2]
    h.a_dofs[2] .-= dt*h.e_dofs[3]
    
end

  
"""
    operatorHE(h, dt)

Push H_E: Equations to be solved
```math
\\begin{aligned}
\\partial_t f + E_1 \\partial_{v_1} f + E_2 \\partial_{v_2} f = 0 &->& V_{new} = V_{old} + dt * E \\\\
\\partial_t E_1 = 0 &->& E_{1,new} = E_{1,old} \\\\
\\partial_t E_2 = 0 &->& E_{2,new} = E_{2,old} \\\\
\\partial_t B + \\partial_{x_1} E_2 = 0 &->& B_{new} = B_{old} - dt \\partial_{x_1} E_2
\\end{aligned}
```
"""
function operatorHE(h :: HamiltonianSplittingSpin, dt :: Float64)

    HH = 0.00022980575
    n_cells = h.kernel_smoother_0.n_dofs

    fill!.(h.part, 0.0)

    hat_v = zeros(Float64, 3, 3)

    S  = zeros(Float64,3)
    St = zeros(Float64,3)
    V  = zeros(Float64,3)
    aa = zeros(Float64,n_cells)
    
    @inbounds for i_part=1:h.particle_group.n_particles

        xi = h.particle_group.particle_array[1, i_part]
        vi = h.particle_group.particle_array[2, i_part]
         
        # Evaluate efields at particle position
        fill!.(h.j_dofs, 0.0)

        add_charge!( h.j_dofs[2], h.kernel_smoother_1, xi, 1.0)
        compute_rderivatives_from_basis!(h.j_dofs[1], h.maxwell_solver, h.j_dofs[2])

        Y  = h.a_dofs[1]'h.j_dofs[1]
        Z  = h.a_dofs[2]'h.j_dofs[1]
        V .= [0, Z, -Y]

        hat_v[1,2] =   Y
        hat_v[1,3] =   Z
        hat_v[2,1] = - Y
        hat_v[3,1] = - Z

        s1 = h.particle_group.particle_array[4, i_part]
        s2 = h.particle_group.particle_array[5, i_part]
        s3 = h.particle_group.particle_array[6, i_part]

        vnorm = norm(V)

        if vnorm > 1e-14
            S .= ( [s1, s2, s3] .+ (sin(dt*vnorm)/vnorm*hat_v 
                   + 0.5* (sin(dt/2*vnorm)/(vnorm/2))^2*hat_v^2)*[s1, s2, s3])
        else
            S .= [s1, s2, s3]
        end

        h.particle_group.particle_array[4, i_part] = S[1]
        h.particle_group.particle_array[5, i_part] = S[2]
        h.particle_group.particle_array[6, i_part] = S[3]

        if vnorm > 1e-14
            St .= (dt .* [s1, s2, s3] + ( 2*(sin(dt*vnorm/2)/vnorm)^2*hat_v 
                                        + 2.0/(vnorm^2)*(dt/2-sin(dt*vnorm)/2/vnorm)*hat_v^2)*[s1, s2, s3])
        else
            St .= dt .* [s1,s2,s3]
        end
        
        wi = get_charge(h.particle_group, i_part) 

        # define part1 and part2
        h.part[1] .+= wi[1]*St[3]*h.j_dofs[2]
        h.part[2] .+= wi[1]*St[2]*h.j_dofs[2]
        
        # update velocity
        fill!(h.j_dofs[2], 0.0)
        add_charge!( h.j_dofs[2], h.kernel_smoother_2, xi, 1.0)
        compute_rderivatives_from_basis!(h.j_dofs[1], h.maxwell_solver, h.j_dofs[2])
        compute_rderivatives_from_basis!(aa, h.maxwell_solver, h.j_dofs[1])
        h.j_dofs[1] .= aa
        vi = vi - HH  * (h.a_dofs[2]'*h.j_dofs[1] * St[2] + h.a_dofs[1]'*h.j_dofs[1] * St[3])

        set_v(h.particle_group, i_part, vi)

    end
    
    # Update bfield
    compute_rderivatives_from_basis!(h.j_dofs[1], h.maxwell_solver,  h.part[1])
    compute_rderivatives_from_basis!(h.j_dofs[2], h.maxwell_solver, -h.part[2])

    h.j_dofs[1] .*= HH
    h.j_dofs[2] .*= HH
    
    compute_e_from_j!( h.e_dofs[2], h.maxwell_solver, h.j_dofs[1], 2) 
    compute_e_from_j!( h.e_dofs[3], h.maxwell_solver, h.j_dofs[2], 2)
        
end

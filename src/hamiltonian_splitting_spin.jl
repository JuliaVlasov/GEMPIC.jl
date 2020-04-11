export SpinHamiltonianSplitting


"""
    HamiltonianSplitting( maxwell_solver,
                          kernel_smoother_0, kernel_smoother_1,
                          particle_group, e_dofs, b_dofs, domain) 

Hamiltonian splitting type for Vlasov-Maxwell

- Integral over the spline function on each interval (order p+1)
- Integral over the spline function on each interval (order p)
- `e_dofs` describing the two components of the electric field
- `b_dofs` describing the magnetic field
- `j_dofs` for kernel representation of current density. 
"""
struct SpinHamiltonianSplitting

    maxwell_solver    :: AbstractMaxwellSolver
    kernel_smoother_0 :: SpinParticleMeshCoupling
    kernel_smoother_1 :: SpinParticleMeshCoupling
    kernel_smoother_2 :: SpinParticleMeshCoupling
    particle_group    :: SpinParticleGroup

    spline_degree     :: Int64
    Lx                :: Float64
    x_min             :: Float64
    delta_x           :: Float64

    cell_integrals_0  :: SVector
    cell_integrals_1  :: SVector

    e_dofs :: Array{Array{Float64,1}}
    a_dofs :: Array{Array{Float64,1}}
    j_dofs :: Array{Array{Float64,1}}
    part1  :: Array{Float64,1}
    part2  :: Array{Float64,1}
    part3  :: Array{Float64,1}
    part4  :: Array{Float64,1}

    function SpinHamiltonianSplitting( maxwell_solver,
                                   kernel_smoother_0,
                                   kernel_smoother_1,
                                   kernel_smoother_2,
                                   particle_group,
                                   e_dofs,
                                   a_dofs,
                                   domain :: Vector{Float64}, nx) 

        # Check that n_dofs is the same for both kernel smoothers.
        @assert kernel_smoother_0.n_dofs == kernel_smoother_1.n_dofs

        j_dofs = [zeros(Float64,kernel_smoother_0.n_dofs) for i in 1:2]
        part1 = zeros(Float64, nx)
        part2 = zeros(Float64, nx)
        part3 = zeros(Float64, nx)
        part4 = zeros(Float64, nx)
        x_min = domain[1]
        Lx    = domain[3]
        spline_degree = 3
        delta_x = Lx/kernel_smoother_1.n_dofs
        
        cell_integrals_1 = SVector{3}([0.5, 2.0, 0.5] ./ 3.0)
        cell_integrals_0 = SVector{4}([1.0,11.0,11.0,1.0] ./ 24.0)

        new( maxwell_solver, kernel_smoother_0, kernel_smoother_1, kernel_smoother_2, 
             particle_group, spline_degree, Lx, x_min, delta_x,
             cell_integrals_0, cell_integrals_1, e_dofs, a_dofs, j_dofs,
            part1, part2, part3, part4)

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
function strang_splitting!( h           :: SpinHamiltonianSplitting,
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
function operatorHp1(h :: SpinHamiltonianSplitting, dt :: Float64)

    n_cells = h.kernel_smoother_0.n_dofs

    fill!(h.j_dofs[1], 0.0)
    fill!(h.j_dofs[2], 0.0)

    for i_part = 1:h.particle_group.n_particles  

       # Read out particle position and velocity
       x_old = get_x(h.particle_group, i_part)
       vi    = get_v(h.particle_group, i_part)
        
       # Then update particle position:  X_new = X_old + dt * V
       x_new = x_old .+ dt * vi

       # Get charge for accumulation of j
       wi     = get_charge(h.particle_group, i_part)
       qoverm = h.particle_group.q_over_m

       add_current_update_v!( h.j_dofs[1], 
                              h.kernel_smoother_1,
                              x_old, 
                              x_new, 
                              wi[1],
                              qoverm, 
                              vi)

       # Accumulate rho for Poisson diagnostics
       #add_charge!( h.j_dofs[2],
        #            h.kernel_smoother_0, 
        #            x_new, 
        #            wi[1])
      
       x_new[1] = mod(x_new[1], h.Lx)
       set_x(h.particle_group, i_part, x_new)
       

    end

    # Update the electric field.
    compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.j_dofs[1], 1)

end

function operatorHp2(h :: SpinHamiltonianSplitting, dt :: Float64)
    
    n_cells = h.kernel_smoother_0.n_dofs
    fill!(h.part1, 0.0)
    fill!(h.part2, 0.0)
    fill!(h.part3, 0.0)
    fill!(h.part4, 0.0)
    

    qm = h.particle_group.q_over_m
    # Update v_1
    for i_part=1:h.particle_group.n_particles
        fill!(h.j_dofs[1], 0.0)
        fill!(h.j_dofs[2], 0.0)
       # Evaluate b at particle position (splines of order p)
       xi    = get_x(h.particle_group, i_part)       
       vi    = get_v(h.particle_group, i_part)
        
        wi = get_charge(h.particle_group, i_part) 
        add_charge!( h.j_dofs[2], h.kernel_smoother_0, xi, 1.0)# R0 
        add_charge!( h.j_dofs[1], h.kernel_smoother_1, xi, 1.0)# R1注意修改最后一个子系统中这个地方，我搞错了B样条的关系
       # values of the derivatives of basis function
        aa = zeros(Float64,n_cells)
        compute_derivatives_from_basis!(aa, h.maxwell_solver, h.j_dofs[1])
        h.j_dofs[1] .= aa
        vi = vi - dt/2*(h.a_dofs[1]'*h.j_dofs[1] * (h.j_dofs[2]'*h.a_dofs[1]))
        vi = vi - dt/2*(h.a_dofs[1]'*h.j_dofs[2] * (h.j_dofs[1]'*h.a_dofs[1]))
        vi = vi - dt/2*(h.a_dofs[2]'*h.j_dofs[1] * (h.j_dofs[2]'*h.a_dofs[2]))
        vi = vi - dt/2*(h.a_dofs[2]'*h.j_dofs[2] * (h.j_dofs[1]'*h.a_dofs[2]))
        
        set_v(h.particle_group, i_part, vi)
        
        # below we solve electric field
        # first define part1 and part2 to be 0 vector
        h.part1 .= h.part1 .+ dt*wi*(h.j_dofs[2]'*h.a_dofs[1])*h.j_dofs[2]
        h.part2 .= h.part2 .+ dt*wi*(h.j_dofs[2]'*h.a_dofs[2])*h.j_dofs[2]
        
        
        
       # Scale vi by weight to combine both factors 
       #for accumulation of integral over j

       

       

    end

    # Update the electric field. Also, we still need to scale with 1/Lx 
    
    

    compute_e_from_j!( h.e_dofs[2], h.maxwell_solver, -h.part1, 2)
    compute_e_from_j!( h.e_dofs[3], h.maxwell_solver, -h.part2, 2)
    # define part3 and part4
    compute_derivatives_from_basis2!(h.part3, h.maxwell_solver, h.a_dofs[1])
    compute_derivatives_from_basis2!(h.part4, h.maxwell_solver, h.a_dofs[2])
    
    compute_e_from_b!( h.e_dofs[2], h.maxwell_solver, dt, h.part3)
    compute_e_from_b!( h.e_dofs[3], h.maxwell_solver, dt, h.part4)
    
    
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
function operatorHB(h :: SpinHamiltonianSplitting, dt :: Float64)
    for i_part=1:h.particle_group.n_particles

       # Evaluate b at particle position (splines of order p)
       
       
       vi    = get_v(h.particle_group, i_part)
        
        xi = get_x(h.particle_group, i_part)
        e1 = evaluate(h.kernel_smoother_1, xi[1], h.e_dofs[1])
        vi = vi + dt  * e1;  
        
       set_v(h.particle_group, i_part, vi)
    end

    h.a_dofs[1] = h.a_dofs[1] .- dt*h.e_dofs[2];
    h.a_dofs[2] = h.a_dofs[2] .- dt*h.e_dofs[3];
    
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
function operatorHE(h :: SpinHamiltonianSplitting, dt :: Float64)
    HH = 0.00022980575
    n_cells = h.kernel_smoother_0.n_dofs
    fill!(h.part1, 0.0)
    fill!(h.part2, 0.0)
    fill!(h.part3, 0.0)
    fill!(h.part4, 0.0)
    hat_v = zeros(Float64, 3, 3)
    
    for i_part=1:h.particle_group.n_particles

       v_new = get_v( h.particle_group, i_part)
        
       # Evaluate efields at particle position
       xi = get_x(h.particle_group, i_part)
        fill!(h.j_dofs[1], 0.0)
        fill!(h.j_dofs[2], 0.0)
        add_charge!( h.j_dofs[2], h.kernel_smoother_1, xi, 1.0)
        compute_derivatives_from_basis!(h.j_dofs[1], h.maxwell_solver, h.j_dofs[2])
        Y = h.a_dofs[1]'*h.j_dofs[1]
        Z = h.a_dofs[2]'*h.j_dofs[1]
        V = [0,Z,-Y];
        hat_v[1,2] = Y;
        hat_v[1,3] = Z;
        hat_v[2,1] = -Y;
        hat_v[3,1] = -Z;
        s1 = get_s1(h.particle_group, i_part)
        s2 = get_s2(h.particle_group, i_part)
        s3 = get_s3(h.particle_group, i_part)
        S = zeros(Float64,3)
        if norm(V)>10^(-14)
            S .= [s1, s2, s3] + ( sin(dt*norm(V))/norm(V)*hat_v + 0.5* (sin(dt/2*norm(V))/(norm(V)/2))^2*hat_v^2 )*[s1, s2, s3]
        else
            S .= [s1, s2, s3]
        end
        set_s1(h.particle_group, i_part, S[1])
        set_s2(h.particle_group, i_part, S[2])
        set_s3(h.particle_group, i_part, S[3])
        St = zeros(Float64,3)
        if norm(V)>10^(-14)
            St .= dt*[s1, s2, s3] + ( 2*(sin(dt*norm(V)/2)/norm(V))^2*hat_v + 2.0/(norm(V)^2)*(dt/2-  sin(dt*norm(V))/2/norm(V)  )*hat_v^2 )*[s1, s2, s3]
        else
            St .= dt*[s1,s2,s3]
        end
        
        wi = get_charge(h.particle_group, i_part) 
        # define part1 and part2
        h.part1 .= h.part1 .+ wi[1]*St[3]*h.j_dofs[2]
        h.part2 .= h.part2 .+ wi[1]*St[2]*h.j_dofs[2]
        #add_charge!( h.part1, h.kernel_smoother_1, xi, wi[1]*St[3])
        #add_charge!( h.part2, h.kernel_smoother_1, xi, wi[1]*St[2])
        
        # update velocity
        fill!(h.j_dofs[2], 0.0)
        add_charge!( h.j_dofs[2], h.kernel_smoother_2, xi, 1.0)
        compute_derivatives_from_basis!(h.j_dofs[1], h.maxwell_solver, h.j_dofs[2])
        aa = zeros(Float64,n_cells)
        compute_derivatives_from_basis!(aa, h.maxwell_solver, h.j_dofs[1])
        h.j_dofs[1] .= aa
        vi = v_new - HH  * (h.a_dofs[2]'*h.j_dofs[1] * St[2] + h.a_dofs[1]'*h.j_dofs[1] * St[3])

        set_v(h.particle_group, i_part, vi)
        

    end
    
    # Update bfield
    aa = zeros(Float64,n_cells)
    bb = zeros(Float64,n_cells)
    compute_derivatives_from_basis!(aa, h.maxwell_solver, h.part1)
    compute_derivatives_from_basis!(bb, h.maxwell_solver, -h.part2)
    
    compute_e_from_j!( h.e_dofs[2], h.maxwell_solver, HH .* aa, 2) # 注意compute_e_from_j里面最后是负号，我们要做相应调整
    compute_e_from_j!( h.e_dofs[3], h.maxwell_solver, HH .* bb, 2)
        
end
  


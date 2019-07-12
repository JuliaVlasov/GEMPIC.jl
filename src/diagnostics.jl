import GEMPIC: get_charge, get_x, add_charge!, get_mass
import GEMPIC: inner_product, evaluate

export solve_poisson!

"""
   solve_poisson!( efield, particle_group, kernel_smoother, maxwell_solver, rho )

Accumulate rho and solve Poisson
 - particle_group : Particles
 - maxwell_solver : Maxwell solver (FEM 1D)
 - kernel_smoother_0 : Particle-Mesh method
 - rho : Charge density
 - efield_dofs : Electric field (1D)
"""
function solve_poisson!( efield_dofs       :: Vector{Float64},
                         particle_group    :: ParticleGroup, 
                         kernel_smoother_0 :: ParticleMeshCoupling, 
                         maxwell_solver    :: Maxwell1DFEM, 
                         rho               :: Vector{Float64}) 

    
    fill!(rho, 0.0)

    for i_part = 1:particle_group.n_particles
       xi = get_x(particle_group, i_part)
       wi = get_charge(particle_group, i_part)
       add_charge!(rho, kernel_smoother_0, xi, wi)
    end

    # Solve Poisson problem
    compute_e_from_rho!( efield_dofs, maxwell_solver, rho )

end


"""
  pic_diagnostics_transfer( particle_group, kernel_smoother_0, 
                            kernel_smoother_1, efield_dofs, transfer)

  Compute ``\\sum_{particles} w_p ( v_1,p e_1(x_p) + v_2,p e_2(x_p)) ``

- particle_group   
- kernel_smoother_0  : Kernel smoother (order p+1)
- kernel_smoother_1  : Kernel smoother (order p)   
- efield_dofs : coefficients of efield

"""
function pic_diagnostics_transfer( particle_group, 
                                   kernel_smoother_0, kernel_smoother_1, 
                                   efield_1_dofs, efield_2_dofs )

    transfer = 0.0
    for i_part = 1:particle_group.n_particles

       xi = get_x( particle_group,i_part )
       wi = get_charge(particle_group, i_part )
       vi = get_v( particle_group, i_part )

       efield_1 = evaluate( kernel_smoother_1, xi[1], efield_1_dofs )
       efield_2 = evaluate( kernel_smoother_0, xi[1], efield_2_dofs )

       transfer += (vi[1] * efield_1 + vi[2] * efield_2) * wi
       
    end

    transfer
    
end

"""
    pic_diagnostics_vvb( particle_group, kernel_smoother_1, bfield_dofs )

Compute ``\\sum_{particles} ( w_p v_1, p b(x_p) v_2, p )``

- particle_group    : particle group object
- kernel_smoother_1 : Kernel smoother (order p)  
- bfield_dofs : coefficients of bfield

"""
function pic_diagnostics_vvb( particle_group, kernel_smoother_1, bfield_dofs )

    vvb = 0.0
    for i_part = 1:particle_group.n_particles

       xi = get_x( particle_group, i_part )
       wi = get_charge( particle_group, i_part )
       vi = get_v( particle_group, i_part )

       bfield = evaluate(kernel_smoother_1, xi[1], bfield_dofs )

       vvb += wi * vi[1] * vi[2] * bfield
     
    end

    vvb

end

"""
    pic_diagnostics_poynting( maxwell_solver, degree, efield_dofs, bfield_dofs, 
                              rho )

Compute ``e^T M_0^{-1}  R^T b``

- maxwell_solver : maxwell solver object
- degree : degree of finite element
- efield_dofs : coefficients of efield
- bfield_dofs : coefficients of bfield

"""
function pic_diagnostics_poynting( maxwell_solver, degree, efield_dofs, 
                                   bfield_dofs )

    scratch = similar(bfield_dofs)
    # Multiply B by M_0^{-1}  R^T
    compute_e_from_b!(scratch, maxwell_solver, 1.0, bfield_dofs )

    inner_product( maxwell_solver, efield_dofs, scratch, degree )

end

  
  
export time_history_diagnostics


"""
    time_history_diagnostics( particle_group, maxwell_solver,
       kernel_smoother_0, kernel_smoother_1,
       time, degree, efield_dofs, bfield_dofs,
       efield_dofs_n, efield_poisson)

Diagnostics for PIC
- particle_group : Particles data
- maxwell_solver : Maxwell solver
- kernel_smoother_0 : Mesh coupling operator
- kernel_smoother_1 : Mesh coupling operator
- time : Time
- efield_dofs : Electric field
- efield_dofs_n : Electric field at half step
- efield_poisson : Electric field from Poisson equation
- rho : charge density
- bfield_dofs : Magnetic field
- degree : Spline degree
"""
function time_history_diagnostics(
       particle_group, maxwell_solver,
       kernel_smoother_0, kernel_smoother_1,
       time, degree, efield_1_dofs, efield_2_dofs, 
       bfield_dofs, efield_1_dofs_n, efield_2_dofs_n, 
       efield_poisson)


    diagnostics = zeros(Float64, 3)

    for i_part=1:particle_group.n_particles
       vi = get_v(   particle_group, i_part)
       wi = get_mass(particle_group, i_part)
       # Kinetic energy
       diagnostics[1] += (vi[1]^2+vi[2]^2)*wi[1]
       # Momentum 1
       diagnostics[2] += vi[1] * wi[1]
       # Momentum 2
       diagnostics[3] += vi[2] * wi[1]
    end

    transfer = pic_diagnostics_transfer( particle_group, kernel_smoother_0, 
                              kernel_smoother_1, efield_1_dofs, efield_2_dofs )

    vvb = pic_diagnostics_vvb( particle_group, kernel_smoother_1, bfield_dofs)

    poynting = pic_diagnostics_poynting( maxwell_solver, degree, 
                                         efield_2_dofs, bfield_dofs )
  
    potential_energy = zeros(Float64, 3)

    potential_energy[1] = inner_product(maxwell_solver, efield_1_dofs, 
                                        efield_1_dofs_n, degree-1 )

    potential_energy[2] = inner_product(maxwell_solver, efield_2_dofs, 
                                        efield_2_dofs_n, degree )

    potential_energy[3] = l2norm_squared( maxwell_solver, bfield_dofs, degree-1 )

    println( """ 
time = $time,  
potential_energy = $potential_energy, 
diagnostics[1] = $(diagnostics[1]),
diagnostics + sum(potential_energy) = $(diagnostics[1] + sum(potential_energy))
diagnostics[2:3] = $(diagnostics[2:3]), 
-transfer+vvb+poynting =  $(-transfer+vvb+poynting),
maximum(abs.(efield_1_dofs .- efield_poisson))) = $(maximum(abs.(efield_1_dofs .- efield_poisson)))
""")

end 

    
export evaluate

"""
    evaluate( kernel_smoother, field_dofs,  xi, n_dofs )

Evaluate the field at points xi

- field_dofs : field value on dofs
- xi : positions where the field is evaluated
"""
function evaluate( kernel_smoother :: ParticleMeshCoupling, 
                   field_dofs :: AbstractArray,  x :: AbstractArray )

    field_grid = similar(x)
    for j in eachindex(x)
        field_grid[j] = evaluate(kernel_smoother, x[j], field_dofs)
    end
    field_grid

end 

"""
    pic_diagnostics_hpi( particle_group,  index, kinetic )

compute v(index)-part of kinetic energy

- particle_group 
- index : velocity component
- kinetic : value of index part of kinetic energy

"""
function pic_diagnostics_hpi( particle_group,  index, kinetic )
    
    kinetic = 0.0
    for i_part = 1:particle_group.n_particles
       vi = get_v(   particle_group, i_part)
       wi = get_mass(particle_group, i_part)
       kinetic += kinetic_local + (vi[index]^2) * wi[1]
    end

    kinetic
    
end 

"""
    eval_derivative_spline( position, xmin, delta_x, n_grid, 
                            field_dofs, degree, derivative )

Compute the spline coefficient of the derivative of some given spline expansion

- position : particle position
- xmin : lower boundary of the domain
- delta_x : step 
- n_grid : number of grid points
- field_dofs : coefficients of spline representation of the field
- degree : degree of spline
- derivative : value of the derivative
"""
function eval_derivative_spline( position, xmin, delta_x, 
                                 n_grid, field_dofs, degree )
    
    der_degree = degree-1
    
    xi = (position[1] - xmin)/delta_x
    index = ceil(xi)
    xi = xi - (index-1)
    index = index - der_degree

    spline_val = uniform_bsplines_eval_basis( der_degree, xi )
    
    derivative = 0.0

    for i1 = 1:degree
       ind = mod(index+i1-2, n_grid)+1
       derivative += spline_val[i1]*(field_dofs[ind]-field_dofs[mod(ind-2, n_grid)+1])
    end

    derivative/delta_x
    
end 

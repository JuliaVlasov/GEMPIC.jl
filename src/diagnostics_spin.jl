using DataFrames

export solve_poisson!

"""
   solve_poisson!( efield, particle_group, kernel_smoother, maxwell_solver, rho )

Accumulate rho and solve Poisson
 - `particle_group` : Particles
 - `maxwell_solver` : Maxwell solver (FEM 1D)
 - `kernel_smoother_0` : Particle-Mesh method
 - `rho` : preallocated array for Charge density
 - `efield_dofs` : spline coefficients of electric field (1D)
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

    compute_e_from_rho!( efield_dofs, maxwell_solver, rho )

end


"""
    pic_diagnostics_transfer( particle_group, kernel_smoother_0, 
                            kernel_smoother_1, efield_dofs, transfer)

Compute ``\\sum_{particles} w_p ( v_1,p e_1(x_p) + v_2,p e_2(x_p)) ``

- `particle_group`   
- `kernel_smoother_0`  : Kernel smoother (order p+1)
- `kernel_smoother_1`  : Kernel smoother (order p)   
- `efield_dofs` : coefficients of efield

"""
function pic_diagnostics_transfer( particle_group, 
                                   kernel_smoother_0, kernel_smoother_1, 
                                   efield_dofs )

    transfer = 0.0
    for i_part = 1:particle_group.n_particles

       xi = get_x( particle_group,i_part )
       wi = get_charge(particle_group, i_part )
       vi = get_v( particle_group, i_part )

       efield_1 = evaluate( kernel_smoother_1, xi[1], efield_dofs[1] )
       efield_2 = evaluate( kernel_smoother_0, xi[1], efield_dofs[2] )

       transfer += (vi[1] * efield_1) * wi
       
    end

    transfer
    
end

"""
    pic_diagnostics_vvb( particle_group, kernel_smoother_1, bfield_dofs )

Compute ``\\sum_{particles} ( w_p v_1, p b(x_p) v_2, p )``

- `particle_group`    : particle group object
- `kernel_smoother_1` : Kernel smoother (order p)  
- `bfield_dofs` : coefficients of bfield

"""
#=
function pic_diagnostics_vvb( particle_group, kernel_smoother_1, afield_dofs )

    vvb = 0.0
    for i_part = 1:particle_group.n_particles

       xi = get_x( particle_group, i_part )
       wi = get_charge( particle_group, i_part )
       vi = get_v( particle_group, i_part )

        bfield = evaluate(kernel_smoother_1, xi[1], afield_dofs[1] )

       vvb += wi * vi[1] *  bfield
     
    end

    vvb

end
=#
"""
    pic_diagnostics_poynting( maxwell_solver, degree, efield_dofs, bfield_dofs, 
                              rho )

Compute ``e^T M_0^{-1}  R^T b``

- `maxwell_solver` : maxwell solver object
- `degree` : degree of finite element
- `efield_dofs` : coefficients of `efield`
- `bfield_dofs` : coefficients of `bfield`

"""
#=
function pic_diagnostics_poynting( maxwell_solver, degree, efield_dofs, 
                                   afield_dofs )

    scratch = similar(bfield_dofs)
    # Multiply B by M_0^{-1}  R^T
    compute_e_from_b!(scratch, maxwell_solver, 1.0, afield_dofs[1] )

    inner_product( maxwell_solver, efield_dofs, scratch, degree )

end
=#
  
export TimeHistoryDiagnostics

"""
    TimeHistoryDiagnostics( particle_group, maxwell_solver, 
                            kernel_smoother_0, kernel_smoother_1 )

Context to save and plot diagnostics

- `particle_group` : Particles data
- `maxwell_solver` : Maxwell solver
- `kernel_smoother_0` : Mesh coupling operator
- `kernel_smoother_1` : Mesh coupling operator
- `data` : DataFrame containing time history values
"""
mutable struct TimeHistoryDiagnostics

    particle_group    :: ParticleGroup
    maxwell_solver    :: Maxwell1DFEM
    kernel_smoother_0 :: ParticleMeshCoupling
    kernel_smoother_1 :: ParticleMeshCoupling
    data              :: DataFrame

    function TimeHistoryDiagnostics( particle_group    :: ParticleGroup,
                                     maxwell_solver    :: Maxwell1DFEM,
                                     kernel_smoother_0 :: ParticleMeshCoupling,
                                     kernel_smoother_1 :: ParticleMeshCoupling)


        data = DataFrame(Time = Float64[],
                         KineticEnergy = Float64[],
                         Kineticspin = Float64[],
                         Momentum1 = Float64[],
                         Momentum2 = Float64[],
                         Momentum3 = Float64[],
                         Momentum4 = Float64[],
                         Momentum5 = Float64[],
                         Momentum6 = Float64[],
                         Momentum7 = Float64[],
                         PotentialEnergyE1 = Float64[],
                         PotentialEnergyE2 = Float64[],
                         PotentialEnergyE3 = Float64[],
                         PotentialEnergyB2 = Float64[],
                         PotentialEnergyB3 = Float64[],
                         Transfer = Float64[],
                         #VVB = Float64[],
                         #Poynting = Float64[],
                         ErrorPoisson = Float64[])


        new( particle_group, maxwell_solver, kernel_smoother_0,
             kernel_smoother_1, data )
    end
end

export write_step!

"""
    write_step!( thdiag, time, degree, efield_dofs, bfield_dofs,
                 efield_dofs_n, efield_poisson)

write diagnostics for PIC
- `time` : Time
- `efield_dofs` : Electric field
- `efield_dofs_n` : Electric field at half step
- `efield_poisson` : Electric field compute from Poisson equation
- `bfield_dofs` : Magnetic field
- `degree` : Spline degree
"""
function write_step!( thdiag :: TimeHistoryDiagnostics,
                      time, degree, efield_dofs,
                      afield_dofs, efield_dofs_n, 
                      efield_poisson, propagator)

    diagnostics = zeros(Float64, 9)
    potential_energy = zeros(Float64, 5)
    HH = 0.00022980575
    for i_part=1:thdiag.particle_group.n_particles
        fill!(propagator.j_dofs[1], 0.0)
        fill!(propagator.j_dofs[2], 0.0)
       xi = get_x(   thdiag.particle_group, i_part)
       vi = get_v(   thdiag.particle_group, i_part)
       s1 = get_s1(   thdiag.particle_group, i_part)
       s2 = get_s2(   thdiag.particle_group, i_part)
       s3 = get_s3(   thdiag.particle_group, i_part) 
       wi = get_mass(thdiag.particle_group, i_part)
       # Kinetic energy
       v2 = evaluate(thdiag.kernel_smoother_0, xi[1], afield_dofs[1])
       v3 = evaluate(thdiag.kernel_smoother_0, xi[1], afield_dofs[2])
       diagnostics[1] += 0.5*(vi[1]^2 + v2[1]^2 + v3[1]^2) * wi[1] 
       add_charge!( propagator.j_dofs[2], propagator.kernel_smoother_1, xi, 1.0)
       compute_derivatives_from_basis!(propagator.j_dofs[1], propagator.maxwell_solver, propagator.j_dofs[2])
       diagnostics[2] += HH * (afield_dofs[2]'*propagator.j_dofs[1]*wi*s2 -  afield_dofs[1]'*propagator.j_dofs[1]*wi*s3)
       # Momentum 1
       diagnostics[3] += xi[1] * wi[1]
       diagnostics[4] += xi[1] * wi[1] * s1
       diagnostics[5] += xi[1] * wi[1] * s2
       diagnostics[6] += xi[1] * wi[1] * s3
       diagnostics[7] += wi[1] * s1
       diagnostics[8] += wi[1] * s2
       diagnostics[9] += wi[1] * s3

    end

    transfer = pic_diagnostics_transfer( thdiag.particle_group, 
        thdiag.kernel_smoother_0, thdiag.kernel_smoother_1, 
        efield_dofs )
#=
    vvb = pic_diagnostics_vvb( thdiag.particle_group, 
    thdiag.kernel_smoother_1, afield_dofs[1])

    poynting = pic_diagnostics_poynting( thdiag.maxwell_solver, degree, 
    efield_dofs[2], afield_dofs[1] )
=#  
    potential_energy[1] = 0.5*inner_product(thdiag.maxwell_solver, 
        efield_dofs[1], efield_dofs[1], degree-1 )
    
    potential_energy[2] = 0.5*inner_product( thdiag.maxwell_solver, 
        efield_dofs[2], efield_dofs[2], degree )
   
    potential_energy[3] = 0.5*inner_product( thdiag.maxwell_solver, 
        efield_dofs[3], efield_dofs[3], degree )
    nn =  thdiag.kernel_smoother_0.n_dofs
    aa = zeros(Float64,nn)
    compute_derivatives_from_basis2!(aa, thdiag.maxwell_solver, afield_dofs[1])
    potential_energy[4] = 0.5*l2norm_squared( thdiag.maxwell_solver, 
        aa, degree-1 )
    bb = zeros(Float64,nn)
    compute_derivatives_from_basis2!(bb, thdiag.maxwell_solver, afield_dofs[2])
    potential_energy[5] = 0.5*l2norm_squared( thdiag.maxwell_solver, 
        bb, degree-1 )

#=
    println( """ 
time = $time,  
potential_energy = $potential_energy, 
diagnostics[1] = $(diagnostics[1]),
diagnostics + sum(potential_energy) = $(diagnostics[1] + sum(potential_energy))
diagnostics[2:3] = $(diagnostics[2:3]), 
-transfer+vvb+poynting =  $(-transfer+vvb+poynting),
maximum(abs.(efield_1_dofs .- efield_poisson))) = $(maximum(abs.(efield_1_dofs .- efield_poisson)))
""")
=#
    push!(thdiag.data, ( time,  
                         diagnostics...,
                         potential_energy..., 
                         transfer,
                         maximum(abs.(efield_dofs[1] .- efield_poisson))))
end 

    
export evaluate

"""
    evaluate( kernel_smoother, field_dofs,  xi, n_dofs )

Evaluate the field at points xi

- `field_dofs` : field value on dofs
- `xi` : positions where the field is evaluated
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

- `particle_group` 
- `index` : velocity component
- `kinetic` : value of index part of kinetic energy

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

- `position` : particle position
- `xmin` : lower boundary of the domain
- `delta_x` : step 
- `n_grid` : number of grid points
- `field_dofs` : coefficients of spline representation of the field
- `degree` : degree of spline
- `derivative` : value of the derivative
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


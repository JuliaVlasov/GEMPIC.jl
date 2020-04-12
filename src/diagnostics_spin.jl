using DataFrames

"""
    pic_diagnostics_transfer( particle_group, kernel_smoother, efield_dofs)

Compute ``\\sum_{particles} w ( v, p e(x_p) ``

- `particle_group`   
- `kernel_smoother` : Kernel smoother
- `efield_dofs` : coefficients of efield

"""
function pic_diagnostics_transfer( particle_group :: ParticleGroup, 
                                   kernel_smoother, efield_dofs :: Vector{Float64} )

    transfer = 0.0

    @inbounds for i_part = 1:particle_group.n_particles

       xi = get_x( particle_group,i_part )
       wi = get_charge(particle_group, i_part )
       vi = get_v( particle_group, i_part )

       efield = evaluate( kernel_smoother, xi, efield_dofs )

       transfer += (vi * efield) * wi
       
    end

    transfer
    
end

export TimeHistoryDiagnosticsSpin

"""
    TimeHistoryDiagnosticsSpin( particle_group, maxwell_solver, 
                            kernel_smoother_0, kernel_smoother_1 )

Context to save and plot diagnostics

- `particle_group` : Particles data
- `maxwell_solver` : Maxwell solver
- `kernel_smoother_0` : Mesh coupling operator
- `kernel_smoother_1` : Mesh coupling operator
- `data` : DataFrame containing time history values
"""
mutable struct TimeHistoryDiagnosticsSpin

    particle_group    :: ParticleGroup
    maxwell_solver    :: Maxwell1DFEM
    kernel_smoother_0 :: ParticleMeshCoupling
    kernel_smoother_1 :: ParticleMeshCoupling
    data              :: DataFrame

    function TimeHistoryDiagnosticsSpin( particle_group :: ParticleGroup,
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
function write_step!( thdiag :: TimeHistoryDiagnosticsSpin,
                      time, degree, efield_dofs,
                      afield_dofs, efield_dofs_n, 
                      efield_poisson, propagator)

    diagnostics = zeros(Float64, 9)
    potential_energy = zeros(Float64, 5)
    HH = 0.00022980575

    nn =  thdiag.kernel_smoother_0.n_dofs
    aa = zeros(Float64,nn)
    bb = zeros(Float64,nn)

    for i_part=1:thdiag.particle_group.n_particles
        fill!(propagator.j_dofs[1], 0.0)
        fill!(propagator.j_dofs[2], 0.0)
        xi = get_x( thdiag.particle_group, i_part)
        vi = get_v( thdiag.particle_group, i_part)
        s1 = get_spin( thdiag.particle_group, i_part, 1)
        s2 = get_spin( thdiag.particle_group, i_part, 2)
        s3 = get_spin( thdiag.particle_group, i_part, 3) 
        wi = get_mass(thdiag.particle_group, i_part)
        # Kinetic energy
        v2 = evaluate(thdiag.kernel_smoother_0, xi[1], afield_dofs[1])
        v3 = evaluate(thdiag.kernel_smoother_0, xi[1], afield_dofs[2])
        diagnostics[1] += 0.5*(vi[1]^2 + v2[1]^2 + v3[1]^2) * wi[1] 
        add_charge!( propagator.j_dofs[2], propagator.kernel_smoother_1, xi, 1.0)
        compute_rderivatives_from_basis!(propagator.j_dofs[1], propagator.maxwell_solver, propagator.j_dofs[2])
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

    transfer = pic_diagnostics_transfer( thdiag.particle_group, thdiag.kernel_smoother_0, efield_dofs[1] )

    potential_energy[1] = 0.5*inner_product( thdiag.maxwell_solver, efield_dofs[1], efield_dofs[1], degree-1 )
    potential_energy[2] = 0.5*inner_product( thdiag.maxwell_solver, efield_dofs[2], efield_dofs[2], degree )
    potential_energy[3] = 0.5*inner_product( thdiag.maxwell_solver, efield_dofs[3], efield_dofs[3], degree )

    compute_lderivatives_from_basis!(aa, thdiag.maxwell_solver, afield_dofs[1])

    potential_energy[4] = 0.5*l2norm_squared( thdiag.maxwell_solver, aa, degree-1 )

    compute_lderivatives_from_basis!(bb, thdiag.maxwell_solver, afield_dofs[2])

    potential_energy[5] = 0.5*l2norm_squared( thdiag.maxwell_solver, bb, degree-1 )

    push!(thdiag.data, ( time,  
                         diagnostics...,
                         potential_energy..., 
                         transfer,
                         maximum(abs.(efield_dofs[1] .- efield_poisson))))
end 

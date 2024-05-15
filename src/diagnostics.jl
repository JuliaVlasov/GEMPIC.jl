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
function solve_poisson!(
    efield_dofs::Vector{Float64},
    particle_group::AbstractParticleGroup,
    kernel_smoother_0::AbstractParticleMeshCoupling,
    maxwell_solver::Maxwell1DFEM,
    rho::Vector{Float64},
)
    fill!(rho, 0.0)

    for i_part in 1:(particle_group.n_particles)
        xi = particle_group.array[1, i_part]
        wi = get_charge(particle_group, i_part)
        add_charge!(rho, kernel_smoother_0, xi, wi)
    end

    return compute_e_from_rho!(efield_dofs, maxwell_solver, rho)
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
function pic_diagnostics_transfer(
    particle_group::ParticleGroup{1,2}, kernel_smoother_0, kernel_smoother_1, efield_dofs
)
    transfer = 0.0::Float64

    for i_part in 1:(particle_group.n_particles)
        xi = particle_group.array[1, i_part]
        wi = get_charge(particle_group, i_part)
        v1 = particle_group.array[2, i_part]
        v2 = particle_group.array[3, i_part]

        efield_1 = evaluate(kernel_smoother_1, xi, efield_dofs[1])
        efield_2 = evaluate(kernel_smoother_0, xi, efield_dofs[2])

        transfer += (v1 * efield_1 + v2 * efield_2) * wi
    end

    return transfer
end

"""
    pic_diagnostics_vvb( particle_group, kernel_smoother_1, bfield_dofs )

Compute ``\\sum_{particles} ( w_p v_1, p b(x_p) v_2, p )``

- `particle_group`    : particle group object
- `kernel_smoother_1` : Kernel smoother (order p)  
- `bfield_dofs` : coefficients of bfield

"""
function pic_diagnostics_vvb(particle_group, kernel_smoother_1, bfield_dofs)
    vvb = 0.0::Float64

    for i_part in 1:(particle_group.n_particles)
        xi = particle_group.array[1, i_part]
        v1 = particle_group.array[2, i_part]
        v2 = particle_group.array[3, i_part]
        wi = particle_group.array[4, i_part]
        wi *= particle_group.charge
        wi *= particle_group.common_weight

        bfield = evaluate(kernel_smoother_1, xi[1], bfield_dofs)

        vvb += wi * v1 * v2 * bfield
    end

    return vvb
end

"""
    pic_diagnostics_poynting( maxwell_solver, degree, efield_dofs, bfield_dofs, 
                              rho )

Compute ``e^T M_0^{-1}  R^T b``

- `maxwell_solver` : maxwell solver object
- `degree` : degree of finite element
- `efield_dofs` : coefficients of `efield`
- `bfield_dofs` : coefficients of `bfield`

"""
function pic_diagnostics_poynting(maxwell_solver, degree, efield_dofs, bfield_dofs)
    scratch = similar(bfield_dofs)
    # Multiply B by M_0^{-1}  R^T
    compute_e_from_b!(scratch, maxwell_solver, 1.0, bfield_dofs)

    return inner_product(maxwell_solver, efield_dofs, scratch, degree)
end

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
struct TimeHistoryDiagnostics
    particle_group::ParticleGroup
    maxwell_solver::Maxwell1DFEM
    kernel_smoother_0::ParticleMeshCoupling1D
    kernel_smoother_1::ParticleMeshCoupling1D
    data::DataFrame
    diagnostics::Vector{Float64}
    potential_energy::Vector{Float64}

    function TimeHistoryDiagnostics(
            particle_group::ParticleGroup{D, V},
        maxwell_solver::Maxwell1DFEM,
        kernel_smoother_0::ParticleMeshCoupling1D,
        kernel_smoother_1::ParticleMeshCoupling1D,
        ) where {D, V}
        data = DataFrame(;
            Time=Float64[],
            KineticEnergy=Float64[],
            Momentum1=Float64[],
            Momentum2=Float64[],
            PotentialEnergyE1=Float64[],
            PotentialEnergyE2=Float64[],
            PotentialEnergyB3=Float64[],
            Transfer=Float64[],
            VVB=Float64[],
            Poynting=Float64[],
            ErrorPoisson=Float64[],
        )

        diagnostics = zeros(Float64, 3)
        potential_energy = zeros(Float64, 3)

        return new(
            particle_group,
            maxwell_solver,
            kernel_smoother_0,
            kernel_smoother_1,
            data,
            diagnostics,
            potential_energy,
        )
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
function write_step!(
    thdiag::TimeHistoryDiagnostics,
    time,
    degree,
    efield_dofs,
    bfield_dofs,
    efield_dofs_n,
    efield_poisson,
)
    fill!(thdiag.diagnostics, 0.0)
    fill!(thdiag.potential_energy, 0.0)

    for i_part in 1:(thdiag.particle_group.n_particles)
        v1 = thdiag.particle_group.array[2, i_part]
        v2 = thdiag.particle_group.array[3, i_part]
        wi = thdiag.particle_group.array[4, i_part]
        wi *= thdiag.particle_group.mass
        wi *= thdiag.particle_group.common_weight

        # Kinetic energy
        thdiag.diagnostics[1] += (v1^2 + v2^2) * wi
        # Momentum 1
        thdiag.diagnostics[2] += v1 * wi
        # Momentum 2
        thdiag.diagnostics[3] += v2 * wi
    end

    transfer = pic_diagnostics_transfer(
        thdiag.particle_group,
        thdiag.kernel_smoother_0,
        thdiag.kernel_smoother_1,
        efield_dofs,
    )

    vvb = pic_diagnostics_vvb(thdiag.particle_group, thdiag.kernel_smoother_1, bfield_dofs)

    poynting = pic_diagnostics_poynting(
        thdiag.maxwell_solver, degree, efield_dofs[2], bfield_dofs
    )

    thdiag.potential_energy[1] = inner_product(
        thdiag.maxwell_solver, efield_dofs[1], efield_dofs_n[1], degree - 1
    )

    thdiag.potential_energy[2] = inner_product(
        thdiag.maxwell_solver, efield_dofs[2], efield_dofs_n[2], degree
    )

    thdiag.potential_energy[3] = l2norm_squared(
        thdiag.maxwell_solver, bfield_dofs, degree - 1
    )

    return push!(
        thdiag.data,
        (
            time,
            thdiag.diagnostics...,
            thdiag.potential_energy...,
            transfer,
            vvb,
            poynting,
            maximum(abs.(efield_dofs[1] .- efield_poisson)),
        ),
    )
end

"""
    evaluate( kernel_smoother, field_dofs,  xi, n_dofs )

Evaluate the field at points xi

- `field_dofs` : field value on dofs
- `xi` : positions where the field is evaluated
"""
function evaluate(
    kernel_smoother::ParticleMeshCoupling1D, field_dofs::AbstractArray, x::AbstractArray
)
    field_grid = similar(x)
    for j in eachindex(x)
        field_grid[j] = evaluate(kernel_smoother, x[j], field_dofs)
    end
    return field_grid
end

"""
    pic_diagnostics_hpi( particle_group,  index, kinetic )

compute v(index)-part of kinetic energy

- `particle_group` 
- `index` : velocity component
- `kinetic` : value of index part of kinetic energy

"""
function pic_diagnostics_hpi(particle_group, index, kinetic)
    kinetic = 0.0
    for i_part in 1:(particle_group.n_particles)
        vi = get_v(particle_group, i_part)
        wi = get_mass(particle_group, i_part)
        kinetic += kinetic_local + (vi[index]^2) * wi[1]
    end

    return kinetic
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
function eval_derivative_spline(position, xmin, delta_x, n_grid, field_dofs, degree)
    der_degree = degree - 1

    xi = (position[1] - xmin) / delta_x
    index = ceil(xi)
    xi = xi - (index - 1)
    index = index - der_degree

    spline_val = uniform_bsplines_eval_basis(der_degree, xi)

    derivative = 0.0

    for i1 in 1:degree
        ind = mod(index + i1 - 2, n_grid) + 1
        derivative +=
            spline_val[i1] * (field_dofs[ind] - field_dofs[mod(ind - 2, n_grid) + 1])
    end

    return derivative / delta_x
end

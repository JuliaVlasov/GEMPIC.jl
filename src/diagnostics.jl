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

function solve_poisson!(
    efield_dofs::Vector{Float64},
    particle_group::ParticleGroup{1,V},
    kernel_smoother_0::ParticleMeshCoupling1D,
    maxwell_solver::Maxwell1DFEM,
    rho::Vector{Float64},
) where {V}
    deposit_charge!(rho, particle_group, kernel_smoother_0)
    return compute_e_from_rho!(efield_dofs, maxwell_solver, rho)
end

function solve_poisson!(
    efield_dofs::Vector{Vector{Float64}},
    particle_group::ParticleGroup{2,2},
    kernel::ParticleMeshCoupling2D,
    maxwell_solver::TwoDMaxwell,
    rho::Vector{Float64},
)
    deposit_charge!(rho, particle_group, kernel)
    return compute_e_from_rho!(efield_dofs, maxwell_solver.poisson, rho)
end

function deposit_charge!(
    rho::Vector{Float64},
    particle_group::ParticleGroup{2,2},
    kernel::ParticleMeshCoupling2D,
)
    fill!(rho, 0.0)
    for i in 1:particle_group.n_particles
        add_charge!(
            rho,
            kernel,
            particle_group.array[1, i],
            particle_group.array[2, i],
            get_charge(particle_group, i),
        )
    end
    return rho
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
    pic_diagnostics_transfer(particle_group::ParticleGroup{1,1}, ...)

Compute the longitudinal particle-to-electric-field transfer for 1d1v.
"""
function pic_diagnostics_transfer(
    particle_group::ParticleGroup{1,1}, kernel_smoother_0, kernel_smoother_1, efield_dofs
)
    transfer = 0.0
    for i_part in 1:particle_group.n_particles
        xi = particle_group.array[1, i_part]
        wi = particle_group.array[3, i_part] * particle_group.charge
        wi *= particle_group.common_weight
        v1 = particle_group.array[2, i_part]
        efield_1 = evaluate(kernel_smoother_1, xi, efield_dofs[1])
        transfer += v1 * efield_1 * wi
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
function pic_diagnostics_vvb(
    particle_group::ParticleGroup{1,2}, kernel_smoother_1, bfield_dofs
)
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

pic_diagnostics_vvb(::ParticleGroup{1,1}, kernel_smoother_1, bfield_dofs) = 0.0

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

export TimeHistoryDiagnostics2D

struct TimeHistoryDiagnostics2D
    particle_group::ParticleGroup{2,2}
    maxwell_solver::TwoDMaxwell
    kernel_ex::ParticleMeshCoupling2DAnisotropic
    kernel_ey::ParticleMeshCoupling2DAnisotropic
    density_kernel::ParticleMeshCoupling2D
    rho::Vector{Float64}
    gauss_rho::Vector{Float64}
    data::DataFrame

    function TimeHistoryDiagnostics2D(
        particle_group::ParticleGroup{2,2}, maxwell_solver::TwoDMaxwell
    )
        degree = maxwell_solver.s_deg_0
        grid = maxwell_solver.mesh
        kernel_ex =
            ParticleMeshCoupling2DAnisotropic(grid, degree - 1, degree, :galerkin)
        kernel_ey =
            ParticleMeshCoupling2DAnisotropic(grid, degree, degree - 1, :galerkin)
        density_kernel = ParticleMeshCoupling2D(
            particle_group, grid, degree, :galerkin
        )
        n_dofs = grid.nx * grid.ny
        data = DataFrame(
            Time=Float64[],
            KineticEnergy=Float64[],
            MomentumX=Float64[],
            MomentumY=Float64[],
            ElectricEnergyEx=Float64[],
            ElectricEnergyEy=Float64[],
            MagneticEnergyBz=Float64[],
            Transfer=Float64[],
            GaussLawError=Float64[],
        )
        return new(
            particle_group,
            maxwell_solver,
            kernel_ex,
            kernel_ey,
            density_kernel,
            zeros(n_dofs),
            zeros(n_dofs),
            data,
        )
    end
end

function write_step!(
    diagnostics::TimeHistoryDiagnostics2D,
    time,
    efield_dofs::Vector{Vector{Float64}},
    bfield_dofs::Vector{Vector{Float64}},
)
    particles = diagnostics.particle_group
    kinetic = momentum_x = momentum_y = transfer = 0.0
    for i in 1:particles.n_particles
        vx, vy = particles.array[3, i], particles.array[4, i]
        mass = get_mass(particles, i)
        charge = get_charge(particles, i)
        kinetic += (vx^2 + vy^2) * mass
        momentum_x += vx * mass
        momentum_y += vy * mass
        x, y = particles.array[1, i], particles.array[2, i]
        transfer += charge * (
            vx * evaluate(diagnostics.kernel_ex, x, y, efield_dofs[1]) +
            vy * evaluate(diagnostics.kernel_ey, x, y, efield_dofs[2])
        )
    end

    deposit_charge!(diagnostics.rho, particles, diagnostics.density_kernel)
    compute_rho_from_e!(diagnostics.gauss_rho, diagnostics.maxwell_solver, efield_dofs)
    mean_charge = sum(diagnostics.rho) / length(diagnostics.rho)
    gauss_law_error = maximum(
        abs.(diagnostics.gauss_rho .- (diagnostics.rho .- mean_charge))
    )
    electric_energy_ex = inner_product(
        diagnostics.maxwell_solver, efield_dofs[1], efield_dofs[1], 1, 1
    )
    electric_energy_ey = inner_product(
        diagnostics.maxwell_solver, efield_dofs[2], efield_dofs[2], 2, 1
    )
    magnetic_energy_bz = inner_product(
        diagnostics.maxwell_solver, bfield_dofs[3], bfield_dofs[3], 3, 3
    )

    return push!(
        diagnostics.data,
        (
            time,
            kinetic,
            momentum_x,
            momentum_y,
            electric_energy_ex,
            electric_energy_ey,
            magnetic_energy_bz,
            transfer,
            gauss_law_error,
        ),
    )
end

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
    if !(thdiag.particle_group.dims in ((1, 1), (1, 2)))
        throw(ArgumentError("Diagnostics are implemented for 1d1v and 1d2v particle groups"))
    end
    fill!(thdiag.diagnostics, 0.0)
    fill!(thdiag.potential_energy, 0.0)

    for i_part in 1:(thdiag.particle_group.n_particles)
        v1 = thdiag.particle_group.array[2, i_part]
        v2 = thdiag.particle_group.dims == (1, 2) ?
             thdiag.particle_group.array[3, i_part] : 0.0
        weight_index = thdiag.particle_group.dims[1] + thdiag.particle_group.dims[2] + 1
        wi = thdiag.particle_group.array[weight_index, i_part]
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

    poynting = if thdiag.particle_group.dims == (1, 1)
        0.0
    else
        pic_diagnostics_poynting(
            thdiag.maxwell_solver, degree, efield_dofs[2], bfield_dofs
        )
    end

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
        kinetic += (vi[index]^2) * wi[1]
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

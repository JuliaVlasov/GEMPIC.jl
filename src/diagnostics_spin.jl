using DataFrames

"""
    pic_diagnostics_transfer( particle_group, kernel_smoother, efield_dofs)

Compute ``\\sum_{particles} w  p e(x_p) ``

- `particle_group`   
- `kernel_smoother` : Kernel smoother
- `efield_dofs` : coefficients of efield

"""
function pic_diagnostics_transfer(
    particle_group::ParticleGroup{1,1}, kernel_smoother, efield_dofs::Vector{Float64}
)
    transfer = 0.0

    @inbounds for i_part in 1:(particle_group.n_particles)
        xi = particle_group.array[1, i_part]
        vi = particle_group.array[2, i_part]
        wi =
            particle_group.array[3, i_part] *
            particle_group.charge *
            particle_group.common_weight
        efield = evaluate(kernel_smoother, xi, efield_dofs)

        transfer += (vi * efield) * wi
    end

    return transfer
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
struct TimeHistoryDiagnosticsSpin
    particle_group::ParticleGroup
    maxwell_solver::Maxwell1DFEM
    kernel_smoother_0::ParticleMeshCoupling1D
    kernel_smoother_1::ParticleMeshCoupling1D
    data::DataFrame
    diagnostics::Vector{Float64}
    potential_energy::Vector{Float64}

    function TimeHistoryDiagnosticsSpin(
        particle_group::ParticleGroup,
        maxwell_solver::Maxwell1DFEM,
        kernel_smoother_0::ParticleMeshCoupling1D,
        kernel_smoother_1::ParticleMeshCoupling1D,
    )
        data = DataFrame(;
            Time=Float64[],
            KineticEnergy=Float64[],
            Kineticspin=Float64[],
            Momentum1=Float64[],
            Momentum2=Float64[],
            Momentum3=Float64[],
            Momentum4=Float64[],
            Momentum5=Float64[],
            Momentum6=Float64[],
            Momentum7=Float64[],
            Momentum8=Float64[],
            Momentum9=Float64[],
            Momentum10=Float64[],
            PotentialEnergyE1=Float64[],
            PotentialEnergyE2=Float64[],
            PotentialEnergyE3=Float64[],
            PotentialEnergyB2=Float64[],
            PotentialEnergyB3=Float64[],
            Transfer=Float64[],
            ErrorPoisson=Float64[],
        )

        diagnostics = zeros(Float64, 12)
        potential_energy = zeros(Float64, 5)

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
- `afield_dofs[1]` : Magnetic Potential Ay
- `afield_dofs[2]` : Magnetic Potential Az
- `efield_dofs[1]` : Longitudinal Electric field Ex
- `efield_poisson` : Electric field compute from Poisson equation
- `efield_dofs[2]` : Ey
- `efield_dofs[3]` : Ez
- `degree` : Spline degree
"""
function write_step!(
    thdiag::TimeHistoryDiagnosticsSpin,
    time,
    degree,
    efield_dofs,
    afield_dofs,
    efield_dofs_n,
    efield_poisson,
    propagator,
)
    nn = thdiag.kernel_smoother_0.n_dofs
    tmp = zeros(Float64, nn)

    fill!(thdiag.diagnostics, 0.0)
    fill!(thdiag.potential_energy, 0.0)

    for i_part in 1:(thdiag.particle_group.n_particles)
        fill!(propagator.j_dofs[1], 0.0)
        fill!(propagator.j_dofs[2], 0.0)

        xi = thdiag.particle_group.array[1, i_part]
        vi = thdiag.particle_group.array[2, i_part]
        wi =
            thdiag.particle_group.array[3, i_part] *
            thdiag.particle_group.mass *
            thdiag.particle_group.common_weight
        s1 = thdiag.particle_group.array[4, i_part]
        s2 = thdiag.particle_group.array[5, i_part]
        s3 = thdiag.particle_group.array[6, i_part]
        # Kinetic energy
        v2 = evaluate(thdiag.kernel_smoother_0, xi[1], afield_dofs[1])
        v3 = evaluate(thdiag.kernel_smoother_0, xi[1], afield_dofs[2])
        thdiag.diagnostics[1] += 0.5 * (vi[1]^2 + v2[1]^2 + v3[1]^2) * wi[1]
        add_charge!(propagator.j_dofs[2], propagator.kernel_smoother_1, xi, 1.0)
        compute_rderivatives_from_basis!(
            propagator.j_dofs[1], propagator.maxwell_solver, propagator.j_dofs[2]
        )
        thdiag.diagnostics[2] +=
            propagator.HH * (
                afield_dofs[2]' * propagator.j_dofs[1] * wi * s2 -
                afield_dofs[1]' * propagator.j_dofs[1] * wi * s3
            )

        # \int (x f) dx dp ds 
        thdiag.diagnostics[3] += xi[1] * wi[1]
        # \int (x s f) dx dp ds 
        thdiag.diagnostics[4] += xi[1] * wi[1] * s1
        thdiag.diagnostics[5] += xi[1] * wi[1] * s2
        thdiag.diagnostics[6] += xi[1] * wi[1] * s3
        # \int (s f) dx dp ds
        thdiag.diagnostics[7] += wi[1] * s1
        thdiag.diagnostics[8] += wi[1] * s2
        thdiag.diagnostics[9] += wi[1] * s3
        # \int (s x B) dx dp ds 
        thdiag.diagnostics[10] += (
            afield_dofs[1]' * propagator.j_dofs[1] * wi * s2 +
            afield_dofs[2]' * propagator.j_dofs[1] * wi * s3
        )
        thdiag.diagnostics[11] += afield_dofs[1]' * propagator.j_dofs[1] * wi * s1 * (-1.0)
        thdiag.diagnostics[12] += afield_dofs[2]' * propagator.j_dofs[1] * wi * s1 * (-1.0)
    end

    transfer = pic_diagnostics_transfer(
        thdiag.particle_group, thdiag.kernel_smoother_0, efield_dofs[1]
    )

    thdiag.potential_energy[1] =
        0.5 *
        inner_product(thdiag.maxwell_solver, efield_dofs[1], efield_dofs[1], degree - 1)
    thdiag.potential_energy[2] =
        0.5 * inner_product(thdiag.maxwell_solver, efield_dofs[2], efield_dofs[2], degree)
    thdiag.potential_energy[3] =
        0.5 * inner_product(thdiag.maxwell_solver, efield_dofs[3], efield_dofs[3], degree)

    compute_lderivatives_from_basis!(tmp, thdiag.maxwell_solver, afield_dofs[1])

    thdiag.potential_energy[4] =
        0.5 * l2norm_squared(thdiag.maxwell_solver, tmp, degree - 1)

    compute_lderivatives_from_basis!(tmp, thdiag.maxwell_solver, afield_dofs[2])

    thdiag.potential_energy[5] =
        0.5 * l2norm_squared(thdiag.maxwell_solver, tmp, degree - 1)

    return push!(
        thdiag.data,
        (
            time,
            thdiag.diagnostics...,
            thdiag.potential_energy...,
            transfer,
            maximum(abs.(efield_dofs[1] .- efield_poisson)),
        ),
    )
end

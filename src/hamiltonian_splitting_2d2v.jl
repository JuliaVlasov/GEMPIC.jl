struct ParticleWorkspace2D
    current_x::Vector{Float64}
    current_y::Vector{Float64}
    kernel_x::ParticleMeshCoupling2DAnisotropic
    kernel_y::ParticleMeshCoupling2DAnisotropic
    kernel_b::ParticleMeshCoupling2DAnisotropic
    kernel_e1::ParticleMeshCoupling2DAnisotropic
    kernel_e2::ParticleMeshCoupling2DAnisotropic
end

"""
    HamiltonianSplitting2D(maxwell, particles, efield, bfield)

Hamiltonian splitting for the closed 2d2v electromagnetic subsystem. It evolves
`Ex`, `Ey` and `Bz`; the uncoupled `Ez`, `Bx` and `By` components remain zero
when initialized to zero.
"""
struct HamiltonianSplitting2D <: AbstractSplitting
    maxwell_solver::TwoDMaxwell
    particle_group::ParticleGroup{2,2}
    e_dofs::Vector{Vector{Float64}}
    b_dofs::Vector{Vector{Float64}}
    current_dofs::Vector{Vector{Float64}}
    chunks::Vector{UnitRange{Int}}
    workspaces::Vector{ParticleWorkspace2D}
    quad_x::Vector{Float64}
    quad_w::Vector{Float64}

    function HamiltonianSplitting2D(
        maxwell_solver::TwoDMaxwell,
        particle_group::ParticleGroup{2,2},
        e_dofs::Vector{Vector{Float64}},
        b_dofs::Vector{Vector{Float64}},
    )
        particle_group.n_particles > 0 ||
            throw(ArgumentError("2d2v splitting requires at least one particle"))
        length(e_dofs) == 3 && length(b_dofs) == 3 ||
            throw(ArgumentError("2d Maxwell fields must have three components"))

        nxny = maxwell_solver.mesh.nx * maxwell_solver.mesh.ny
        chunksize = cld(particle_group.n_particles, nthreads())
        chunks = [
            first:min(first + chunksize - 1, particle_group.n_particles) for
            first in 1:chunksize:particle_group.n_particles
        ]
        degree = maxwell_solver.s_deg_0
        grid = maxwell_solver.mesh
        workspaces = [
            ParticleWorkspace2D(
                zeros(nxny),
                zeros(nxny),
                ParticleMeshCoupling2DAnisotropic(
                    grid, degree - 1, degree, :galerkin
                ),
                ParticleMeshCoupling2DAnisotropic(
                    grid, degree, degree - 1, :galerkin
                ),
                ParticleMeshCoupling2DAnisotropic(
                    grid, degree - 1, degree - 1, :galerkin
                ),
                ParticleMeshCoupling2DAnisotropic(
                    grid, degree - 1, degree, :galerkin
                ),
                ParticleMeshCoupling2DAnisotropic(
                    grid, degree, degree - 1, :galerkin
                ),
            ) for _ in chunks
        ]
        quad_x, quad_w = gausslegendre(degree + 1)

        return new(
            maxwell_solver,
            particle_group,
            e_dofs,
            b_dofs,
            [zeros(nxny), zeros(nxny)],
            chunks,
            workspaces,
            quad_x,
            quad_w,
        )
    end
end

export HamiltonianSplitting2D

function strang_splitting!(
    h::HamiltonianSplitting2D, dt::Float64, number_steps::Int
)
    for _ in 1:number_steps
        operatorHB(h, 0.5dt)
        operatorHE(h, 0.5dt)
        operatorHp2(h, 0.5dt)
        operatorHp1(h, dt)
        operatorHp2(h, 0.5dt)
        operatorHE(h, 0.5dt)
        operatorHB(h, 0.5dt)
    end
    return h
end

function operatorHB(h::HamiltonianSplitting2D, dt::Float64)
    return compute_e_from_b!(h.e_dofs, h.maxwell_solver, dt, h.b_dofs)
end

function operatorHE(h::HamiltonianSplitting2D, dt::Float64)
    particles = h.particle_group.array
    qm = h.particle_group.q_over_m
    @sync for (chunk, workspace) in zip(h.chunks, h.workspaces)
        @spawn begin
            for i in chunk
                x, y = particles[1, i], particles[2, i]
                particles[3, i] += dt * qm * evaluate(workspace.kernel_e1, x, y, h.e_dofs[1])
                particles[4, i] += dt * qm * evaluate(workspace.kernel_e2, x, y, h.e_dofs[2])
            end
        end
    end
    return compute_b_from_e!(h.b_dofs, h.maxwell_solver, dt, h.e_dofs)
end

function operatorHp1(h::HamiltonianSplitting2D, dt::Float64)
    fill!(h.current_dofs[1], 0.0)
    particles = h.particle_group.array
    grid = h.maxwell_solver.mesh
    qm = h.particle_group.q_over_m

    @sync for (chunk, workspace) in zip(h.chunks, h.workspaces)
        @spawn begin
            fill!(workspace.current_x, 0.0)
            for i in chunk
                x, y = particles[1, i], particles[2, i]
                x_new = x + dt * particles[3, i]
                qweight = particles[5, i] * h.particle_group.charge *
                          h.particle_group.common_weight
                vy = deposit_path_2d!(
                    workspace.current_x, workspace.kernel_x, workspace.kernel_b,
                    h.b_dofs[3], h.quad_x, h.quad_w, grid, true, x, y, x_new, y,
                    qweight, qm, particles[4, i], -1.0,
                )
                particles[1, i] = mod(x_new - grid.xmin, grid.dimx) + grid.xmin
                particles[4, i] = vy
            end
        end
    end
    for workspace in h.workspaces
        h.current_dofs[1] .+= workspace.current_x
    end
    return compute_e_from_j!(h.e_dofs[1], h.maxwell_solver, h.current_dofs[1], 1)
end

function operatorHp2(h::HamiltonianSplitting2D, dt::Float64)
    fill!(h.current_dofs[2], 0.0)
    particles = h.particle_group.array
    grid = h.maxwell_solver.mesh
    qm = h.particle_group.q_over_m

    @sync for (chunk, workspace) in zip(h.chunks, h.workspaces)
        @spawn begin
            fill!(workspace.current_y, 0.0)
            for i in chunk
                x, y = particles[1, i], particles[2, i]
                y_new = y + dt * particles[4, i]
                qweight = particles[5, i] * h.particle_group.charge *
                          h.particle_group.common_weight
                vx = deposit_path_2d!(
                    workspace.current_y, workspace.kernel_y, workspace.kernel_b,
                    h.b_dofs[3], h.quad_x, h.quad_w, grid, false, x, y, x, y_new,
                    qweight, qm, particles[3, i], 1.0,
                )
                particles[2, i] = mod(y_new - grid.ymin, grid.dimy) + grid.ymin
                particles[3, i] = vx
            end
        end
    end
    for workspace in h.workspaces
        h.current_dofs[2] .+= workspace.current_y
    end
    return compute_e_from_j!(h.e_dofs[2], h.maxwell_solver, h.current_dofs[2], 2)
end

function deposit_path_2d!(
    current, current_kernel::ParticleMeshCoupling2DAnisotropic,
    magnetic_kernel::ParticleMeshCoupling2DAnisotropic, bfield,
    quad_x, quad_w, grid::TwoDGrid, axis_x::Bool, x0, y0, x1, y1,
    qweight, q_over_m, velocity_transverse, magnetic_sign,
)
    spacing = axis_x ? grid.dx : grid.dy
    lower = axis_x ? grid.xmin : grid.ymin
    start = axis_x ? x0 : y0
    stop = axis_x ? x1 : y1
    delta = (stop - start) / spacing
    delta == 0.0 && return velocity_transverse

    direction = sign(delta)
    u = (start - lower) / spacing
    finish = u + delta
    while direction > 0 ? u < finish : u > finish
        boundary = direction > 0 ? floor(Int, u) + 1 : ceil(Int, u) - 1
        next_u = direction > 0 ? min(finish, boundary) : max(finish, boundary)
        midpoint = 0.5 * (u + next_u)
        halfwidth = 0.5 * (next_u - u)

        @inbounds for k in eachindex(quad_x)
            path_u = midpoint + halfwidth * quad_x[k]
            path_position = lower + spacing * path_u
            xp = axis_x ? path_position : x0
            yp = axis_x ? y0 : path_position

            ix, iy = shape_factor!(current_kernel, xp, yp)
            bx, by = shape_factor!(magnetic_kernel, xp, yp)
            ds = spacing * halfwidth * quad_w[k]
            magnetic_value = 0.0
            nx = grid.nx
            for a in 1:(current_kernel.degree_x + 1),
                b in 1:(current_kernel.degree_y + 1)
                index = mod(ix + a - 2, nx) +
                        mod(iy + b - 2, grid.ny) * nx + 1
                current[index] +=
                    qweight * ds * current_kernel.values_x[a] *
                    current_kernel.values_y[b] * current_kernel.scaling
            end
            for a in 1:(magnetic_kernel.degree_x + 1),
                b in 1:(magnetic_kernel.degree_y + 1)
                index = mod(bx + a - 2, nx) +
                        mod(by + b - 2, grid.ny) * nx + 1
                magnetic_value +=
                    bfield[index] * magnetic_kernel.values_x[a] *
                    magnetic_kernel.values_y[b]
            end
            velocity_transverse += magnetic_sign * q_over_m * magnetic_value * ds
        end
        u = next_u
    end
    return velocity_transverse
end

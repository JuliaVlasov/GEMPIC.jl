using StaticArrays

abstract type AbstractSplitting end

export HamiltonianSplitting

"""
    HamiltonianSplitting( maxwell_solver,
                          kernel_smoother_0, kernel_smoother_1,
                          particle_group, e_dofs, b_dofs) 

Hamiltonian splitting type for Vlasov-Maxwell

- Integral over the spline function on each interval (order p+1)
- Integral over the spline function on each interval (order p)
- `e_dofs` describing the two components of the electric field
- `b_dofs` describing the magnetic field
- `j_dofs` for kernel representation of current density. 
"""
struct HamiltonianSplitting{D,V}
    dims::Tuple{Int64,Int64}
    maxwell_solver::AbstractMaxwellSolver
    kernel_smoother_0::ParticleMeshCoupling1D
    kernel_smoother_1::ParticleMeshCoupling1D
    particle_group::ParticleGroup

    spline_degree::Int
    Lx::Float64
    x_min::Float64
    delta_x::Float64

    cell_integrals_0::SVector
    cell_integrals_1::SVector

    e_dofs::Array{Array{Float64,1}}
    b_dofs::Array{Float64,1}
    j_dofs::Array{Array{Float64,1}}

    chunks::Iterators.PartitionIterator

    function HamiltonianSplitting{D, V}(
        maxwell_solver, kernel_smoother_0, kernel_smoother_1, particle_group, e_dofs, b_dofs
    ) where {D, V}

        dims = (D, V)

        @assert dims == particle_group.dims
        # Check that n_dofs is the same for both kernel smoothers.
        @assert kernel_smoother_0.n_dofs == kernel_smoother_1.n_dofs

        j_dofs = [zeros(Float64, kernel_smoother_0.n_dofs) for i in 1:2]

        x_min = maxwell_solver.xmin
        Lx = maxwell_solver.Lx
        spline_degree = 3
        delta_x = Lx / kernel_smoother_1.n_dofs

        cell_integrals_1 = SVector{3}([0.5, 2.0, 0.5] ./ 3.0)
        cell_integrals_0 = SVector{4}([1.0, 11.0, 11.0, 1.0] ./ 24.0)

        np = particle_group.n_particles
        n_jobs = nthreads()

        @assert np % n_jobs == 0

        chunks = Iterators.partition(1:np, np ÷ n_jobs)

        return new(
            dims,
            maxwell_solver,
            kernel_smoother_0,
            kernel_smoother_1,
            particle_group,
            spline_degree,
            Lx,
            x_min,
            delta_x,
            cell_integrals_0,
            cell_integrals_1,
            e_dofs,
            b_dofs,
            j_dofs,
            chunks,
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
function strang_splitting!(h::HamiltonianSplitting, dt::Float64, number_steps::Int)
    for i_step in 1:number_steps
        operatorHB(h, 0.5dt)
        operatorHE(h, 0.5dt)
        operatorHp2(h, 0.5dt)
        operatorHp1(h, 1.0dt)
        operatorHp2(h, 0.5dt)
        operatorHE(h, 0.5dt)
        operatorHB(h, 0.5dt)
    end
end

export lie_splitting!
"""
    lie_splitting( h, dt, number_steps)

Lie splitting
"""
function lie_splitting!(h::HamiltonianSplitting, dt::Float64, number_steps::Int)
    for i_step in 1:number_steps
        operatorHE(dt)
        operatorHB(dt)
        operatorHp1(dt)
        operatorHp2(dt)
    end
end

export lie_splitting_back!
"""
    lie_splitting_back(h, dt, number_steps)

Lie splitting (oposite ordering)
"""
function lie_splitting_back!(h::HamiltonianSplitting, dt::Float64, number_steps::Int)
    for i_step in 1:number_steps
        operatorHp2(dt)
        operatorHp1(dt)
        operatorHB(dt)
        operatorHE(dt)
    end
end


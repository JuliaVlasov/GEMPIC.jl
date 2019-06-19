module GEMPIC

    include("mesh.jl")
    include("low_level_bsplines.jl")
    include("splinepp.jl")

    # Maxwell solvers
    abstract type AbstractMaxwellSolver end

    include("maxwell_1d_fem.jl")

    # Particle Group
    include("particle_group.jl")

    # Particle-Mesh coupling
    include("particle_mesh_coupling.jl")

    # Particle sampling
    include("particle_sampling.jl")

    # Splittings

    abstract type AbstractSplitting end

    include("hamiltonian_splitting.jl")

end 

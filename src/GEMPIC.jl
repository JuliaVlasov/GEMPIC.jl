module GEMPIC

    # Utilities
    include("mesh.jl")
    include("low_level_bsplines.jl")
    include("splinepp.jl")
    include("distributions.jl")

    # Maxwell solvers
    include("maxwell_1d_fem.jl")

    # Particle Group
    include("particle_group.jl")

    # Particle-Mesh coupling
    include("particle_mesh_coupling.jl")

    # Particle sampling
    include("particle_sampling.jl")
    include("landau_damping.jl")

    # Splittings
    include("hamiltonian_splitting.jl")
    include("hamiltonian_splitting_boris.jl")

    # Diagnostics
    include("diagnostics.jl")

end 

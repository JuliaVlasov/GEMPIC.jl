module GEMPIC

    # Utilities
    include("mesh.jl")
    include("low_level_bsplines.jl")
    include("splinepp.jl")
    include("distributions.jl")
    include("distributions_spin.jl")


    # Maxwell solvers
    include("maxwell_1d_fem.jl")

    # Particle Groups
    include("particle_group.jl")
    include("particle_group_spin.jl")

    # Particle-Mesh coupling
    include("particle_mesh_coupling.jl")
    include("particle_mesh_coupling_spin.jl")

    # Particle sampling
    include("particle_sampling.jl")
    include("particle_sampling_spin.jl")
    include("landau_damping.jl")

    # Splittings
    include("hamiltonian_splitting.jl")
    include("hamiltonian_splitting_boris.jl")
    include("hamiltonian_splitting_spin.jl")

    # Diagnostics
    include("diagnostics.jl")
    include("diagnostics_spin.jl")

end 

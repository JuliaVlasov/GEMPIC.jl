module GEMPIC

import Base.Threads: @sync, @spawn, nthreads, threadid

# Utilities
include("mesh.jl")
include("low_level_bsplines.jl")
include("splinepp.jl")
include("distributions.jl")

# Maxwell solvers
include("maxwell_1d_fem.jl")
include("linear_solver_spline_mass_2d.jl")
include("poisson_2d_fem.jl")
include("maxwell_2d_fem.jl")

# Particle Groups
include("particle_group.jl")

# Particle-Mesh coupling
abstract type AbstractParticleMeshCoupling end
include("particle_mesh_coupling_1d.jl")
include("particle_mesh_coupling_2d.jl")

# Particle sampling
include("particle_sampling.jl")
include("landau_damping.jl")

# Splittings
include("hamiltonian_splitting.jl")
include("hamiltonian_splitting_1d1v.jl")
include("hamiltonian_splitting_1d2v.jl")
include("hamiltonian_splitting_2d3v.jl")
include("hamiltonian_splitting_boris.jl")

# Diagnostics
include("diagnostics.jl")

end

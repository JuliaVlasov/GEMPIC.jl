using Test 
using GEMPIC

#include("test_particle_mesh_coupling_spline_2d.jl")
include("test_vm_1d1v_spin.jl")
include("test_vm_1d2v.jl")
include("test_mesh.jl")
include("test_spline_pp.jl")
include("test_maxwell_1d_fem.jl")
include("test_sampling.jl")
include("test_sampling_spin.jl")
include("test_particle_mesh_coupling_spline_1d.jl")

if Threads.nthreads() == 1 
    include("test_hamiltonian_splitting.jl")
end

include("test_hamiltonian_splitting_spin.jl")
include("test_hamiltonian_splitting_boris.jl")
include("test_save_particles.jl")

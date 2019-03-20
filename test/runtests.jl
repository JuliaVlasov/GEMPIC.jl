using Test

include("test_maxwell_1d_fem.jl")

@testset " Maxwell 1D FEM solver " begin

    mode = 2
    test_maxwell_1d_fem( mode )

end

include("test_particle_mesh_coupling_spline_1d.jl")

@testset " particle-mesh coupling with spline 1d " begin

    test_particle_mesh_coupling_spline_1d( )

end

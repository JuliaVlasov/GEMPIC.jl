using Test

@testset " particle-mesh coupling with spline 1d " begin

    include("test_particle_mesh_coupling_spline_1d.jl")
    test_particle_mesh_coupling_spline_1d( )

end

@testset " Maxwell 1D FEM solver " begin

    include("test_maxwell_1d_fem.jl")
    mode = 2
    test_maxwell_1d_fem( mode )

end


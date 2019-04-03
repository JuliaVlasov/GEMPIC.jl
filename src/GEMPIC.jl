module GEMPIC

    include("mesh.jl")
    include("low_level_bsplines.jl")
    include("splinepp.jl")
    include("maxwell_1d_fem.jl")
    include("particle_group.jl")
    include("particle_mesh_coupling.jl")
    include("particle_sampling.jl")
    include("pic_vm.jl")

end 

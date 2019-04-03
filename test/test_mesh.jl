using GEMPIC

@testset " Mesh " begin


    mesh1d = Mesh( -π, π, 11 )
    mesh2d = Mesh( 0.0, 2π, 11, -6.0, 6.0, 21 )
    mesh3d = Mesh( 0.0, 2π, 11, -π, π, 21, 0.0, 1.0, 32 )

    @test get_x(mesh1d, 1)           ≈ [-π]
    @test get_x(mesh1d, 6)           ≈ [0.0]
    @test get_x(mesh2d, [1,1])       ≈ [0.0, -6.0]
    @test get_x(mesh2d, [6,11])      ≈ [π, 0.0]
    @test get_x(mesh3d, [1,1,1])     ≈ [0.0, -π, 0.0]
    @test get_x(mesh3d, [6, 11, 32]) ≈ [π, 0.0, 1.0]

end

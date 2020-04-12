@testset "VM 1D1V with spin" begin

   σ, μ = 0.02, 0.0
   kx, α = 1.004355, 0.001
   xmin, xmax = 0, 2π/kx
   domain = [xmin, xmax, xmax - xmin]
   nx = 8 
   n_particles = 2
   mesh = Mesh( xmin, xmax, nx)
   spline_degree = 3
   
   df = CosSumGaussianSpin([[kx]], [α], [[σ]], [[μ]] )
   
   particle_group2 = ParticleGroup{1,1}( n_particles, n_spin=3)   
   sampler = ParticleSampler{1,1}( :sobol, n_particles)
   
   sample!(particle_group2, sampler, df, mesh)

   ref =  [3.127970342747129    4.691955514120693; 
           0.01734694403902491 -0.018034876317136343;
           6.249684744808763    6.255940685494258;
           0.0                  0.0; 
           0.0                  0.0; 
           1.0                  1.0] 

   @test particle_group2.particle_array ≈ ref
   
   particle_group = ParticleGroup{1,1}( n_particles, n_spin=3)   
   GEMPIC.set_common_weight(particle_group, (1.0/n_particles))

   for  i_part = 1:n_particles

       x = zeros( 1 )
       v = zeros( 1 )
       s = zeros( 3 )
       w = zeros( 1 )
       x = GEMPIC.get_x(particle_group2, i_part)
       v = GEMPIC.get_v(particle_group2, i_part)

       s1 = GEMPIC.get_spin(particle_group2, i_part, 1)
       s2 = GEMPIC.get_spin(particle_group2, i_part, 2)
       s3 = GEMPIC.get_spin(particle_group2, i_part, 3)

       w = GEMPIC.get_weights(particle_group2, i_part)

       GEMPIC.set_x(particle_group, i_part, x[1])
       GEMPIC.set_v(particle_group, i_part, v[1])
       GEMPIC.set_spin(particle_group, i_part, 1, s1)
       GEMPIC.set_spin(particle_group, i_part, 2, s2)
       GEMPIC.set_spin(particle_group, i_part, 3, s3)
       GEMPIC.set_weights(particle_group, i_part, w[1])

   end
   
   kernel_smoother2 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree-2, :galerkin) 
   kernel_smoother1 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree-1, :galerkin)    
   kernel_smoother0 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree, :galerkin)
   
   rho = zeros(Float64, nx)
   efield_poisson = zeros(Float64, nx)
   
   # Init!ialize the field solver
   maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
   # efield by Poisson
   solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )

   ref = [-0.5811006688338523, -0.8336082614158025, -2.7872165528999013, 1.3710914975949415, -1.3738807953831014, 2.7885897256327588, 0.8341900701471621, 0.5819349851577953]

   @test efield_poisson ≈ ref
   
   # Initialize the arrays for the spline coefficients of the fields
   k0 = 12.0523 
   E0 = 10
   ww = 12.104827940833333

   Ey(x) = E0*cos(k0*x)
   Ez(x) = E0*sin(k0*x)
   Ay(x) = -E0/ww*sin(k0*x)
   Az(x) = E0/ww*cos(k0*x)
   
   efield_dofs = [efield_poisson, zeros(Float64, nx), zeros(Float64, nx)]
   afield_dofs = [zeros(Float64, nx), zeros(Float64, nx)]
   
   l2projection!( efield_dofs[2], maxwell_solver, Ey, spline_degree)
   l2projection!( efield_dofs[3], maxwell_solver, Ez, spline_degree)
   l2projection!( afield_dofs[1], maxwell_solver, Ay, spline_degree)
   l2projection!( afield_dofs[2], maxwell_solver, Az, spline_degree)

   @test efield_dofs[2] ≈ [-0.17074650510231487, 0.17074650468466593, -0.17074650409994616, 0.1707465033481691, -0.17074650242932413, 0.17074650134342606, -0.1707465000904745, 0.1707464986704385]
   @test efield_dofs[3] ≈ [-1.0681800094815476e-5, 1.6022700129154877e-5, -2.136360014868184e-5, 2.670450014853168e-5, -3.204540010485847e-5, 3.738630007371513e-5, -4.272719996293515e-5, 4.806809980886038e-5]
   @test afield_dofs[1] ≈ [8.824412992125462e-7, -1.3236619477267236e-6, 1.7648825950343306e-6, -2.2061032407091485e-6, 2.6473238828043315e-6, -3.088544525898473e-6, 3.5297651624492664e-6, -3.970985795411606e-6]
   @test afield_dofs[2] ≈ [-0.014105653210181833, 0.01410565317567925, -0.014105653127374468, 0.014105653065269032, -0.01410565298936165, 0.014105652899653832, -0.014105652796145382, 0.014105652678833834]

   
   propagator = HamiltonianSplittingSpin( maxwell_solver,
                                      kernel_smoother0, 
                                      kernel_smoother1, 
                                      kernel_smoother2,
                                      particle_group,
                                      efield_dofs,
                                      afield_dofs,
                                      domain, nx);
   
   efield_dofs_n = propagator.e_dofs
   
   thdiag = TimeHistoryDiagnosticsSpin( particle_group, maxwell_solver, 
                           kernel_smoother0, kernel_smoother1 );
   
   steps, Δt = 2, 0.1

   # Strang splitting
   strang_splitting!(propagator, Δt, 1)

   solve_poisson!( efield_poisson, particle_group, 
                   kernel_smoother0, maxwell_solver, rho)

   @test efield_poisson ≈ [-0.5797623802247858, -0.8360479090580495, -2.782777070172533, 1.3790403945715886, -1.381791774035773, 2.7843379708057716, 0.8364431982466378, 0.5805575698671437]
       
   strang_splitting!(propagator, Δt, 1)

   solve_poisson!( efield_poisson, particle_group, 
                   kernel_smoother0, maxwell_solver, rho)

   @test efield_poisson ≈ [-0.5726359838763564, -0.8485012058572868, -2.7596506491700334, 1.4173810569155116, -1.4201086257409283, 2.761425759659844, 0.8486960737209063, 0.5733935743483431]

   @test propagator.particle_group.particle_array ≈ [3.117245017656692 4.702527578518938; 
         -0.1255875828581744 0.12472320092353745; 
          6.249684744808763 6.255940685494258;
         -6.256040292957428e-5 6.112305895910396e-5; 
          3.774570240505298e-9 3.5321392403606012e-9; 
          0.999999998043098 0.9999999981319859] 

   write_step!(thdiag, Δt, spline_degree, 
                       efield_dofs,  afield_dofs,
                       efield_dofs_n, efield_poisson, propagator)


   GEMPIC.save( "particles", 1, particle_group)
       

end

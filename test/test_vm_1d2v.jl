import GEMPIC: get_charge, get_x, add_charge!, get_mass
import GEMPIC: inner_product, evaluate
import GEMPIC: strang_splitting

"""
Simulation of 1d2v Vlasov-Maxwell with simple PIC method, 
periodic boundary conditions, Weibel instability. 
FEM with splines, degree 3 for B and 2 for E
"""

@testset " PIC VM 1D2V " begin

    delta_t         = 0.05
    n_time_steps    = 10
    beta            = 0.0001
    initial_distrib = :cossum_onegaussian
    initial_bfield  = :cos
    
    kx        = 1.25
    alpha     = 0.0
    v_thermal = [0.2,  0.005773502691896]
    v_mean    = [0.0, 0.0]
    
    ng_x   = 32
    x1_min = 0.0
    x1_max = 5.02654824574
    
    n_particles    = 100000
    sampling_case  = :sobol
    symmetric      = true
    splitting_case = :symplectic
    spline_degree  = 3
    
    mesh   = Mesh( x1_min, x1_max, ng_x)
    domain = [x1_min, x1_max, x1_max - x1_min ]

    beta_cos_k(x) = beta * cos(2π * x / domain[3]) 
    beta_sin_k(x) = beta * sin(2π * x / domain[3]) 
    
    n_total_particles = n_particles
    degree_smoother   = spline_degree

    df = CosGaussian( (1,2), 1, 1, hcat([kx]), [alpha], 
                      hcat(v_thermal), hcat(v_mean), 0.0 )
    
    sampler = ParticleSampler( sampling_case, symmetric, (1,2), n_particles)
    
    for propagator in [:symplectic, :boris]
    
        # Initialize the particles   (mass and charge set to 1.0)
        particle_group = ParticleGroup{1,2}( n_particles, 1.0, 1.0, 1)
        
        # Init!ialize the field solver
        maxwell_solver = Maxwell1DFEM(domain, ng_x, spline_degree)
        
        kernel_smoother1 = ParticleMeshCoupling( domain, [ng_x], n_particles, 
                     spline_degree, :galerkin)
        
        kernel_smoother0 = ParticleMeshCoupling( domain, [ng_x], n_particles, 
                     spline_degree, :galerkin)
    
        # Initialize the arrays for the spline coefficients of the fields
        efield1_dofs = zeros(Float64, ng_x)
        efield2_dofs = zeros(Float64, ng_x)
        bfield_dofs  = zeros(Float64, ng_x)
    
        # Initialize the time-splitting propagator
        if splitting_case == :symplectic

            propagator = HamiltonianSplitting( maxwell_solver,
                                              kernel_smoother0, 
                                              kernel_smoother1, 
                                              particle_group,
                                              efield1_dofs, 
                                              efield2_dofs, 
                                              bfield_dofs,
                                              domain[1], 
                                              domain[3]    )

           efield_1_dofs_n = propagator.e_dofs_1
           efield_2_dofs_n = propagator.e_dofs_2

        elseif splitting_case == :boris

            propagator = HamiltonianSplittingBoris( maxwell_solver,
                                                   kernel_smoother0, 
                                                   kernel_smoother1, 
                                                   particle_group,
                                                   efield1_dofs, 
                                                   efield2_dofs, 
                                                   bfield_dofs,
                                                   domain[1], 
                                                   domain[3]    )

           efield_1_dofs_n = propagator.e_dofs_1_mid
           efield_2_dofs_n = propagator.e_dofs_2_mid

        end


        sample( sampler, particle_group, df, mesh )

        # Set the initial fields
        rho = zeros(Float64, ng_x)
        efield_poisson = zeros(Float64, ng_x)

        # efield 1 by Poisson
        solve_poisson!( efield1_dofs, particle_group, kernel_smoother0, maxwell_solver, rho )

        # bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))

        if initial_bfield == :cos
            l2projection!( bfield_dofs, maxwell_solver, beta_cos_k, degree_smoother-1)
        else
            l2projection!( bfield_dofs, maxwell_solver, beta_sin_k, degree_smoother-1)
        end
       
        # In case we use the Boris propagator, we need to initialize 
        # the staggering used in the scheme.
        if splitting_case == :boris
           staggering( propagator, delta_t )
        end

        solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )

        time_history_diagnostics( particle_group, maxwell_solver,
                                  kernel_smoother0, kernel_smoother1, 
                                  0.0,
                                  degree_smoother, 
                                  efield1_dofs, efield2_dofs, bfield_dofs,
                                  efield_1_dofs_n, efield_2_dofs_n, efield_poisson)

        for j = 1:n_time_steps # loop over time

           # Strang splitting
           strang_splitting(propagator, delta_t, 1)

           # Diagnostics
           solve_poisson!( efield_poisson, particle_group, 
                           kernel_smoother0, maxwell_solver, rho)

           time_history_diagnostics( particle_group, maxwell_solver,
                                     kernel_smoother0, kernel_smoother1, 
                                     j * delta_t,
                                     degree_smoother, 
                                     efield1_dofs, efield2_dofs, bfield_dofs,
                                     efield_1_dofs_n, efield_2_dofs_n, efield_poisson)

        end

        # Compute final rho
        fill!(rho, 0.0)
        for i_part = 1:particle_group.n_particles
           xi = get_x(particle_group, i_part)
           wi = get_charge( particle_group, i_part)
           add_charge!(rho, kernel_smoother0, xi, wi)
        end

    end

end


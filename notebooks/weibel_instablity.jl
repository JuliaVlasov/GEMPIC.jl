# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     comment_magics: false
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# # Vlasov–Maxwell in 1D2V
#
# Starting from the Vlasov equation for a particle with charge `q`, and mass `m` given in
#
# $x = (x,0,0)$ and $v = (v_1,v_2,0)$ 
#
# as well as
#
# $$ E(x, t) = (E_1, E_2, 0), \qquad B(x, t) = (0, 0, B3), \qquad f(x, v, t) = f(x, v_1, v_2, t),$$
#
# so that the Vlasov equation takes the form
#
# $$
# \frac{\partial f(x,v,t)}{\partial t} 
# + v_1 \frac{\partial f(x,v,t)}{\partial x} 
# + \frac{q}{m} \Big[ 􏰋\mathbf{E}(x,t)+B_3(x,t)􏰉 \begin{pmatrix}v2\\-v_1\end{pmatrix}􏰊􏰌 \Big] \cdot \nabla f(x,v,t) = 0,
# $$
#
# while Maxwell’s equations become
# $$
# \begin{eqnarray}
# \frac{\partial E_1(x, t)}{\partial t} &=& −J_1(x), \\
# \frac{\partial E_2(x, t)}{\partial t} &=& −\frac{\partial B(x, t)}{\partial x} − J_2(x), \\
# \frac{\partial B(x, t)}{\partial t} &=& −\frac{\partial E_2(x, t)}{\partial x} \\
# \frac{\partial E_1(x,t)}{\partial x} &=& \rho + \rho_B
# \end{eqnarray}
# $$
# with sources given by
# $$
# \rho  = q􏰘 \int f dv, \qquad J_1 =q \int f v_1dv, \qquad  J_2 = q \int f v_2 dv 
# $$
# Note that $\nabla \cdot B = 0$ is manifest.

# # Weibel instability
#
# We study a reduced 1d2v model with a perturbation along $x_1$, a magnetic field along $x_3$ and electric fields along the $x_1$ and $x_2$ directions. Moreover, we assume that the distribution function is independent of $v_3$. The initial distribution and fields are of the form
# $$
# f(x,\mathbf{v},t=0)=\frac{1}{2\pi\sigma_1\sigma_2} \exp \Big( - \frac{1}{2} \big( \frac{v_1^2}{\sigma_1^2} + \frac{v_2^2}{\sigma_2^2} \big) \Big) ( 1 + \alpha cos(kx)), \qquad x \in [0, 2\pi / k),
# $$

# $$ B_3(x, t = 0) = \beta \cos(kx), \qquad E_2(x, t = 0) = 0,$$
#
# and $E_1(x, t = 0)$ is computed from Poisson’s equation.
#
# $$ \sigma_1 = 0.02 / \sqrt{2} $$
# $$ \sigma_2 = \sqrt{12} \sigma_1 $$
#
#
# $ k = 1.25, α = 0 $ and $ \beta = −10^{−4}$. 

using Pkg

Pkg.update();

# +
"""
    pic_vm_1d2v( )
Simulation of 1d2v Vlasov-Maxwell with simple PIC method, 
periodic boundary conditions, Weibel instability. 
FEM with splines, degree 3 for B and 2 for E
"""
function pic_vm_1d2v()

    delta_t         = 0.05
    n_time_steps    = 10
    beta            = 0.0001
    initial_bfield  = :cos
    
    kx        = hcat([1.25])
    alpha     = [0.0]
    v_thermal = hcat([0.2,  0.005773502691896])
    v_mean    = hcat([0.0, 0.0])
    
    ng_x   = 32
    x1_min = 0.0
    x1_max = 2π / kx[1]
    
    n_particles    = 100000
    sampling_case  = :sobol
    symmetric      = true
    splitting_case = :symplectic
    spline_degree  = 3
    
    mesh   = Mesh( x1_min, x1_max, ng_x)
    domain = [x1_min, x1_max, x1_max - x1_min ]
    
    beta_cos_k(x) = beta * cos(2π * x / domain[3]) 
    beta_sin_k(x) = beta * sin(2π * x / domain[3]) 
    
    df = CosSumOneGaussian( (1,2), 1, 1, kx, alpha, v_thermal, v_mean, 1.0 )
    sampler = ParticleSampler( sampling_case, symmetric, (1,2), n_particles)

end
# -

df = pic_vm_1d2v()
x = LinRange(0, 2π/1.25, 32) |> collect
v = LinRange(-6, 6, 64) |> collect
contour(x, v, df( x, v))

# +
    
   
    

    

   
    
    
    
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
# -



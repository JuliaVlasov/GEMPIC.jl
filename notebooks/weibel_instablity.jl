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
#     display_name: Julia 1.1.1
#     language: julia
#     name: julia-1.1
# ---

using ProgressMeter, Plots

include("../src/mesh.jl")
include("../src/distributions.jl")
include("../src/low_level_bsplines.jl")
include("../src/splinepp.jl")
include("../src/maxwell_1d_fem.jl")


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

# +
"""
    pic_vm_1d2v( )
Simulation of 1d2v Vlasov-Maxwell with simple PIC method, 
periodic boundary conditions, Weibel instability. 
FEM with splines, degree 3 for B and 2 for E
"""

Δt    = 0.05
steps = 10
β     = 0.0001

k  = 1.25
α  = 0.0
σ  = [0.2,  0.005773502691896]
μ  = [0.0, 0.0]

nx   = 32
xmin = 0.0
xmax = 2π / k

n_particles    = 100000
sampling_case  = :sobol
symmetric      = true
splitting_case = :symplectic
spline_degree  = 3

mesh   = Mesh( xmin, xmax, nx)
domain = [xmin, xmax, xmax - xmin ]

beta_cos_k(x) = β * cos(2π * x / domain[3]) 
beta_sin_k(x) = beta * sin(2π * x / domain[3]) 

df = CosSumGaussian{1,2}([[k]], [α], [σ], [μ] )
# -

v1min, v1max, nv1 = -0.1, 0.1, 64
v2min, v2max, nv2 = -1.0, 1.0, 64
v1 = LinRange(v1min, v1max, nv1) |> collect
v2 = LinRange(v2min, v2max, nv2) |> collect
f = zeros(Float64,(nv1,nv2))
for i in eachindex(v1), j in eachindex(v2)
    f[i,j] = eval_v_density(df, [v1[i],v2[j]])
end

contour(v1, v2, f)

include("../src/particle_group.jl")


?ParticleGroup

 # Initialize the particles   (mass and charge set to 1.0)
mass, charge = 1.0, 1.0
particle_group = ParticleGroup{1,2}( n_particles, mass, charge, 1)   

include("../src/particle_sampling.jl")

?ParticleSampler

sampler = ParticleSampler{1,2}( sampling_case, symmetric, n_particles)

?sample

sample!(  particle_group, sampler, df, mesh)

?get_x

get_x(particle_group, 1)

get_v(particle_group, 2)

xp = Vector{Float64}[]
for i in 1:n_particles
    push!(xp, vcat(get_x(particle_group,i), 
            get_v(particle_group,i),
            get_weights(particle_group,i)))
end

?scatter

p = vcat(xp[1:100:100000]'...)

scatter(p[:,1], p[:,3], markersize=1)

include("../src/particle_mesh_coupling.jl")

?ParticleMeshCoupling

include("../src/hamiltonian_splitting.jl")

df = pic_vm_1d2v()
x = LinRange(0, 2π/1.25, 32) |> collect
v = LinRange(-6, 6, 64) |> collect
contour(x, v, df( x, v))

# +
    
   
    

    

   
    
    
        
       
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



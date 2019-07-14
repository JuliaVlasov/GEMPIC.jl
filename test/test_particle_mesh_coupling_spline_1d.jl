@testset " particle-mesh coupling with spline 1d " begin

using  GEMPIC
import GEMPIC: set_common_weight, set_x, set_v, set_weights
import GEMPIC: get_x, get_v, get_charge, add_charge!, add_charge_pp!
import GEMPIC: add_current_update_v!, add_current_update_v_pp!
import GEMPIC: b_to_pp, evaluate, evaluate_pp

n_cells       = 10            # Number of cells
n_particles   = 4             # Number of particles
spline_degree = 3             # Spline degree
domain        = [0.0, 2.0]    # x_min, x_max
x_vec = [0.1, 0.65, 0.7, 1.5] # Particle positions
v_vec = [1.5  -3.0  0.0  6.0; 
         0.0   0.5  0.0  0.0]'

particle_group = ParticleGroup{1,2}( n_particles, 1.0, 1.0, 1)

set_common_weight(particle_group, 1.0/n_particles)

for i_part = 1:n_particles
    set_x(particle_group, i_part, x_vec[i_part])
    set_weights(particle_group, i_part, 1.0)
    set_v(particle_group, i_part, v_vec[i_part,:])
end

values_grid = zeros(Float64,(4,1,4))

values_grid[:,1,1] .= [ 2.0833333333333332e-002,  
                        0.47916666666666663, 
                        0.47916666666666663, 
                        2.0833333333333332e-002]

values_grid[:,1,3] .= values_grid[:,1,1]   

values_grid[:,1,4] .= values_grid[:,1,1] 

values_grid[:,1,2] .= [7.0312500000000000e-002,  
                       0.61197916666666663,
                       0.31510416666666663,        
                       2.6041666666666665E-003 ]

kernel = ParticleMeshCoupling( domain, [n_cells], n_particles, 
             spline_degree, :collocation)

# Accumulate rho
rho_dofs  = zeros(Float64, n_cells)
rho_dofs1 = zeros(Float64, n_cells)
for i_part = 1:n_particles
    xi = get_x(particle_group, i_part)
    wi = get_charge(particle_group, i_part)
    add_charge!(rho_dofs, kernel, xi, wi)
    add_charge_pp!(rho_dofs1, kernel, xi, wi)
end

rho_dofs_ref        = zeros(Float64, n_cells)
rho_dofs_ref[8:10] .= values_grid[1:3,1,1]
rho_dofs_ref[1]     = values_grid[4,1,1]
rho_dofs_ref[1:4]  .= rho_dofs_ref[1:4] + values_grid[:,1,2] + values_grid[:,1,3]
rho_dofs_ref[5:8]  .= rho_dofs_ref[5:8] + values_grid[:,1,4]
rho_dofs_ref       .= rho_dofs_ref/n_particles * n_cells/domain[2]

@test maximum(abs.(rho_dofs  .- rho_dofs_ref)) < 1e-15
@test maximum(abs.(rho_dofs1 .- rho_dofs_ref)) < 1e-15

j_dofs  = zeros(Float64, n_cells)
j_dofs1 = zeros(Float64, n_cells)
b_dofs  = zeros(Float64, n_cells)

for i_part = 1:n_particles 

    xi    = get_x(particle_group, i_part)
    wi    = get_charge(particle_group, i_part)
    vi    = get_v(particle_group, i_part)
    vi1   = vi
    x_new = xi .+ vi[1]/10.0

    vi  = add_current_update_v!(    j_dofs,  kernel, xi, x_new, wi, 1.0, b_dofs, vi )
    vi1 = add_current_update_v_pp!( j_dofs1, kernel, xi, x_new, wi, 1.0, b_dofs, vi1 )
     
end

j_dofs_ref = [ 2.4617513020833336e-002,   
               4.0690104166666692e-005, 
                                   0.0,
                                   0.0, 
                                   0.0, 
                                   0.0,  
                                   0.0,  
               6.5104166666666674e-004,
               4.6834309895833329e-002,    
                   0.11535644531249999]

j_dofs_ref .= j_dofs_ref .+ [ -0.16219075520833331, 
                              -0.16219075520833331,
                               -2.52685546875e-002, 
                          -4.0690104166666692e-005, 
                                               0.0, 
                                               0.0, 
                                               0.0, 
                                               0.0, 
                          -4.0690104166666692e-005,  
                               -2.52685546875e-002]

j_dofs_ref .= j_dofs_ref .+ [ 6.5104166666666696e-004,  
                                                  0.0, 
                                                  0.0,
                                                  0.0, 
                              6.5104166666666674e-004, 
                              5.0130208333333329e-002,
                                  0.19986979166666666, 
                                  0.24869791666666663,
                                  0.19986979166666666, 
                              5.0130208333333329e-002]
 
@test maximum(abs.(j_dofs  .- j_dofs_ref)) < 1e-15
@test maximum(abs.(j_dofs1 .- j_dofs_ref)) < 1e-15
  
rho_dofs_pp = b_to_pp(kernel.spline_pp, n_cells, rho_dofs)

particle_values  = zeros(Float64, 4)
particle_values1 = zeros(Float64, 4)
   
for i_part = 1:n_particles
    xi = get_x(particle_group, i_part)
    particle_values[i_part]  = evaluate(    kernel, xi[1], rho_dofs)    
    particle_values1[i_part] = evaluate_pp( kernel, xi[1], rho_dofs_pp)
end

particle_values_ref = [ 1.1560058593749998,
                        2.3149278428819446,
                        2.2656250000000000,
                        1.1512586805555554 ] 

particle_values_ref ./= domain[2]

@test maximum(abs.(particle_values  .- particle_values_ref)) < 1e-15
@test maximum(abs.(particle_values1 .- particle_values_ref)) < 1e-15

end

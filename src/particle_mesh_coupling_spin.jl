"""  
    add_current_update_v_pp!( j_dofs, p, position_old, position_new, 
                              marker_charge)

Add current for one particle and update v 
(according to `H_p1` part in Hamiltonian splitting)
"""
function add_current_update_v_pp!( j_dofs        :: AbstractArray, 
                                   p             :: ParticleMeshCoupling1D, 
                                   position_old  :: Float64, 
                                   position_new  :: Float64, 
                                   marker_charge :: Float64)

    # Read out particle position and velocity
    # Compute index_old, the index of the last DoF on the grid the particle
    # contributes to, and r_old, its position (normalized to cell size one).

    xi = (position_old[1] - p.xmin) / p.delta_x[1]
    index_old = trunc(Int, xi)
    r_old = xi - index_old

    # Compute the new box index index_new and normalized position r_old.
    xi = (position_new[1] - p.xmin) / p.delta_x[1]
    index_new = trunc(Int, xi)
    r_new = xi - index_new
 
    if index_old == index_new
         
        update_jv_pp!(j_dofs, p, r_old, r_new, index_old, marker_charge)

    elseif index_old < index_new

        update_jv_pp!(j_dofs, p, r_old, 1.0, index_old, marker_charge)

        update_jv_pp!(j_dofs, p, 0.0, r_new, index_new, marker_charge)

        for ind = index_old+1:index_new-1
            update_jv_pp!(j_dofs, p, 0.0, 1.0, ind, marker_charge)
        end

    else

        update_jv_pp!(j_dofs, p, 1.0, r_new, index_new, marker_charge)
        update_jv_pp!(j_dofs, p, r_old, 0.0, index_old, marker_charge)

        for ind = index_new+1:index_old-1
             update_jv_pp!(j_dofs, p, 1.0, 0.0, ind, marker_charge)
        end

    end
  

end

"""
    update_jv_pp!( j_dofs, p, lower, upper, index, marker_charge)

Helper function for `add_current_update_v`.
"""
function update_jv_pp!( j_dofs         :: AbstractArray,
                        p              :: ParticleMeshCoupling1D, 
                        lower          :: Float64, 
                        upper          :: Float64, 
                        index          :: Int64, 
                        marker_charge  :: Float64)
                        

   n_cells = p.n_grid[1]

   # Evaluation of the primitive integral at the lower 
   # and upper bound of the gridcell

   horner_primitive_1d(p.spline_val,      p.spline_degree, p.spline_pp.poly_coeffs_fp, lower) 
   horner_primitive_1d(p.spline_val_more, p.spline_degree, p.spline_pp.poly_coeffs_fp, upper) 

   p.spline_val .= (p.spline_val_more .- p.spline_val) .* (p.delta_x[1]) 
   
   ind = 1
   @inbounds for i_grid = index - p.spline_degree:index
       i_mod = mod(i_grid, n_cells ) + 1
       j_dofs[i_mod] = j_dofs[i_mod] + (marker_charge*p.spline_val[ind]* p.scaling)
       ind = ind + 1
   end

end


"""
    add_current_update_v!( j_dofs, p, position_old, position_new, marker_charge) 

Add current for one particle and update v (according to ``H_{p1}``
part in Hamiltonian splitting)

- Read out particle position and velocity
- Compute index_old, the index of the last DoF on the grid the 
particle contributes to, and `r_old`, its position (normalized to cell size one).

"""
function add_current_update_v!( j_dofs        :: AbstractArray,
                                p             :: ParticleMeshCoupling1D, 
                                position_old  :: Float64, 
                                position_new  :: Float64, 
                                marker_charge :: Float64) 


    xi = (position_old - p.xmin) / p.delta_x[1]
    index_old = trunc(Int,xi)
    r_old = xi - index_old

    # Compute the new box index index_new and normalized position r_old.

    xi = (position_new - p.xmin) / p.delta_x[1]
    index_new = trunc(Int, xi)
    r_new = xi - index_new
 
    if index_old == index_new

        if r_old < r_new
            update_jv!(j_dofs, p, r_old, r_new, index_old, marker_charge, 1.0)
        else
            update_jv!(j_dofs, p, r_new, r_old, index_old, marker_charge, -1.0)
        end

    elseif index_old < index_new

        update_jv!(j_dofs, p, r_old, 1.0, index_old, marker_charge, 1.0)

        update_jv!(j_dofs, p, 0.0, r_new, index_new, marker_charge, 1.0)

        for ind = index_old+1:index_new-1
            update_jv!(j_dofs, p, 0.0, 1.0, ind, marker_charge, 1.0)
        end

    else

        update_jv!( j_dofs, p, r_new, 1.0, index_new, marker_charge, -1.0)
        update_jv!( j_dofs, p, 0.0, r_old, index_old, marker_charge, -1.0)

        for ind = index_new+1:index_old-1
            update_jv!( j_dofs, p, 0.0, 1.0, ind, marker_charge, -1.0)
        end

     end    

end

"""
    update_jv!(j_dofs, p, lower, upper, index, marker_charge, qoverm, sign)

Helper function for `add_current_update_v`.
"""
function update_jv!(j_dofs        :: AbstractArray, 
                    p             :: ParticleMeshCoupling1D, 
                    lower         :: Float64, 
                    upper         :: Float64, 
                    index         :: Int64, 
                    marker_charge :: Float64, 
                    sign          :: Float64)
                    

   n_cells = p.n_grid[1]

   c1 = 0.5 * (upper-lower)
   c2 = 0.5 * (upper+lower)

   uniform_bsplines_eval_basis!(p.spline_val, p.spline_degree, c1 * p.quad_x[1]+c2)

   p.spline_val .*= p.quad_w[1] * c1

   @inbounds for j = 2:p.n_quad_points

       uniform_bsplines_eval_basis!( p.spline_val_more, p.spline_degree, c1 * p.quad_x[j]+c2) 

       p.spline_val .+= p.spline_val_more .* p.quad_w[j] .* c1

   end

   p.spline_val .*= sign * p.delta_x[1]
   
   ind = 1
   for i_grid = index - p.spline_degree:index
      i_mod = mod(i_grid, n_cells ) + 1
      j_dofs[i_mod] += marker_charge * p.spline_val[ind] * p.scaling
      ind = ind + 1
   end


end

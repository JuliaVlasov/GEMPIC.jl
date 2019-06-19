export ParticleMeshCoupling

"""
Kernel smoother with splines of arbitrary degree placed on a uniform mesh.
Spline with index i starts at point i

- Value of grid spacing along both directions.
- Definition of the domain: domain(1:2) = x1_min, x1_max
- Number of particles of underlying PIC method (processor local)
- Degree of smoothing kernel spline
- Number of intervals where spline non zero (spline_degree + 1)
- Scaling factor depending on whether Galerkin or collocation
- Number of quadrature points
- scratch data for spline evaluation
- more scratch data for spline evaluation
- quadrature weights and points
"""
mutable struct ParticleMeshCoupling

    dims            :: Int
    domain          :: Vector{Float64}
    delta_x         :: Float64 
    n_grid          :: Vector{Int} 
    n_dofs          :: Int
    no_particles    :: Int
    spline_degree   :: Int
    n_span          :: Int
    scaling         :: Float64
    n_quad_points   :: Int
    spline_val      :: Vector{Float64}
    spline_val_more :: Vector{Float64}
    quad_x          :: Vector{Float64}
    quad_w          :: Vector{Float64}
    spline_pp       :: SplinePP

    function ParticleMeshCoupling( domain         :: Vector{Float64}, 
                                   n_grid         :: Vector{Int64}, 
                                   no_particles   :: Int, 
                                   spline_degree  :: Int, 
                                   smoothing_type :: Symbol )
        dims    = 1
        n_dofs  = prod(n_grid)
        delta_x = (domain[2]-domain[1])/n_grid[1]
        n_span  = spline_degree + 1

        if smoothing_type == :collocation
            scaling = 1.0/delta_x
        elseif smoothing_type == :galerkin
            scaling = 1.0
        else
            throw(ArgumentError( """
                  Smoothing Type $smoothing_type not implemented 
                  for kernel_smoother_spline_1d.
            """))
        end
    
        n_quad_points = (spline_degree+2)÷2

        spline_val      = zeros(Float64, n_span)
        spline_val_more = zeros(Float64, n_span)

        quad_x, quad_w = gausslegendre(n_quad_points)

        spline_pp = SplinePP( spline_degree, n_grid[1])

        new(dims, domain, delta_x, n_grid, n_dofs,
            no_particles, spline_degree,
            n_span, scaling, n_quad_points, spline_val,
            spline_val_more, quad_x, quad_w, spline_pp)

  end
    
     
end
     

"""
Add charge of one particle
- p             : kernel smoother object
- position      : Position of the particle
- marker_charge : Particle weights time charge
- rho_dofs      : Coefficient vector of the charge distribution
"""
function add_charge_pp(p, position, marker_charge, rho_dofs)
    
    xi    = (position[1] - p.domain[1])/p.delta_x
    index = floor(Int64, xi)+1
    xi    = xi    - (index-1)
    index = index - p.spline_degree

    for i in eachindex(p.spline_val)

        p.spline_val[i] = horner_1d(p.spline_degree, 
                                       p.spline_pp.poly_coeffs, xi, i)

    end

    for i = 1:p.n_span
       index1d = mod(index+i-2, p.n_grid[1]) + 1
       rho_dofs[index1d] = rho_dofs[index1d] + (marker_charge * 
                                           p.spline_val[i] * p.scaling)
    end

end


"""
Add charge of one particle
- p              : kernel smoother object
- position       : Position of the particle
- marker_charge  : Particle weights time charge
- rho_dofs       : Coefficient vector of the charge distribution
"""
function add_charge(p             :: ParticleMeshCoupling,
                    position      :: Vector{Float64}, 
                    marker_charge :: Float64, 
                    rho_dofs      :: Vector{Float64})

    xi    = (position[1] - p.domain[1])/p.delta_x[1]
    index = floor(Int64, xi)+1
    xi    = xi - (index-1)
    index = index - p.spline_degree

    p.spline_val .= uniform_bsplines_eval_basis(p.spline_degree, xi)

    for i = 1:p.n_span
       index1d = mod(index+i-2,p.n_grid[1]) + 1
       rho_dofs[index1d] += marker_charge * p.spline_val[i] * p.scaling
    end

end

"""  
Add current for one particle and update v 
(according to H_p1 part in Hamiltonian splitting)
"""
function add_current_update_v_pp(p, position_old, position_new, 
                                 marker_charge, qoverm, bfield_dofs, 
                                 vi, j_dofs)

    # Read out particle position and velocity
    # Compute index_old, the index of the last DoF on the grid the particle
    # contributes to, and r_old, its position (normalized to cell size one).

    xi = (position_old[1] - p.domain[1]) / p.delta_x[1]
    index_old = floor(Int64, xi)
    r_old = xi - index_old

    # Compute the new box index index_new and normalized position r_old.
    xi = (position_new[1] - p.domain[1]) / p.delta_x[1]
    index_new = floor(Int64, xi)
    r_new = xi - index_new
 
    if index_old == index_new
         
        update_jv_pp(p, r_old, r_new, index_old, marker_charge,
                     qoverm,  vi[2], j_dofs, bfield_dofs)

    elseif index_old < index_new

        update_jv_pp(p, r_old, 1.0, index_old, marker_charge,
                     qoverm, vi[2], j_dofs, bfield_dofs)

        update_jv_pp(p, 0.0, r_new, index_new, marker_charge,
                     qoverm, vi[2], j_dofs, bfield_dofs)

        for ind = index_old+1:index_new-1
            update_jv_pp(p, 0.0, 1.0, ind, marker_charge,
                          qoverm, vi[2], j_dofs, bfield_dofs)
        end

    else

        update_jv_pp(p, 1.0, r_new, index_new, marker_charge, qoverm, 
                     vi[2], j_dofs, bfield_dofs)

        update_jv_pp(p, r_old, 0.0, index_old, marker_charge, qoverm,
                     vi[2], j_dofs, bfield_dofs)

        for ind = index_new+1:index_old-1
             update_jv_pp(p, 1.0, 0.0, ind, marker_charge, qoverm,
                          vi[2], j_dofs, bfield_dofs)
        end

    end
  

end

"""
Helper function for \a add_current_update_v.
"""
function update_jv_pp( p :: ParticleMeshCoupling , lower, upper, index, 
                       marker_charge, 
                       qoverm, vi, j_dofs :: Vector{Float64}, 
                       bfield_dofs :: Vector{Float64})

   n_cells = p.n_grid[1]

   # Evaluation of the primitive integral at the lower 
   # and upper bound of the gridcell

   horner_primitive_1d(p.spline_val,      p.spline_degree, p.spline_pp.poly_coeffs_fp, lower) 
   horner_primitive_1d(p.spline_val_more, p.spline_degree, p.spline_pp.poly_coeffs_fp, upper) 

   p.spline_val .= (p.spline_val_more .- p.spline_val) .* (p.delta_x[1]) 
   
   ind = 1
   for i_grid = index - p.spline_degree:index
       i_mod = mod(i_grid, n_cells ) + 1
       j_dofs[i_mod] = j_dofs[i_mod] + (marker_charge*p.spline_val[ind]* p.scaling)
       vi  = vi - qoverm * p.spline_val[ind]*bfield_dofs[i_mod]
       ind = ind + 1
   end

end


"""
Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
"""
function add_current_update_v(p             :: ParticleMeshCoupling, 
                              position_old  :: Vector{Float64}, 
                              position_new  :: Vector{Float64}, 
                              marker_charge :: Float64, 
                              qoverm        :: Float64, 
                              bfield_dofs   :: Vector{Float64}, 
                              vi            :: Vector{Float64}, 
                              j_dofs        :: Vector{Float64})

    # Read out particle position and velocity
    # Compute index_old, the index of the last DoF on the grid the 
    # particle contributes to, and r_old, its position 
    # (normalized to cell size one).

    xi = (position_old[1] - p.domain[1]) / p.delta_x[1]
    index_old = floor(Int64,xi)
    r_old = xi - index_old

    # Compute the new box index index_new and normalized position r_old.
    xi = (position_new[1] - p.domain[1]) / p.delta_x[1]
    index_new = floor(Int64, xi)
    r_new = xi - index_new
 
    if index_old == index_new

        if r_old < r_new
            update_jv(p, r_old, r_new, index_old, marker_charge, 
                      qoverm, 1.0, vi[2], j_dofs, bfield_dofs)
        else
            update_jv(p, r_new, r_old, index_old, marker_charge, qoverm, 
                      -1.0, vi[2], j_dofs, bfield_dofs)
        end

    elseif index_old < index_new

        update_jv(p, r_old, 1.0, index_old, marker_charge, 
                  qoverm, 1.0, vi[2], j_dofs, bfield_dofs)

        update_jv(p, 0.0, r_new, index_new, marker_charge, 
                  qoverm, 1.0, vi[2], j_dofs, bfield_dofs)

        for ind = index_old+1:index_new-1
            update_jv(p, 0.0, 1.0, ind, marker_charge, 
                      qoverm, 1.0, vi[2], j_dofs, bfield_dofs)
        end

    else

        update_jv(p, r_new, 1.0, index_new, marker_charge, qoverm, 
                  -1.0, vi[2], j_dofs, bfield_dofs)
        update_jv(p, 0.0, r_old, index_old, marker_charge, qoverm, 
                  -1.0, vi[2], j_dofs, bfield_dofs)

        for ind = index_new+1:index_old-1
            update_jv(p, 0.0, 1.0, ind, marker_charge, qoverm, 
                      -1.0, vi[2], j_dofs, bfield_dofs)
        end

     end    

end

"""
Helper function for \a add_current_update_v.
"""
function update_jv(p, lower, upper, index, marker_charge, qoverm, 
                   sign, vi, j_dofs, bfield_dofs)

   n_cells = p.n_grid[1]

   c1 = 0.5 * (upper-lower)
   c2 = 0.5 * (upper+lower)

   p.spline_val .= uniform_bsplines_eval_basis(p.spline_degree, 
                                               c1 * p.quad_x[1]+c2)

   p.spline_val .*= p.quad_w[1] * c1


   for j = 2:p.n_quad_points

       p.spline_val_more .= uniform_bsplines_eval_basis(p.spline_degree, 
                                                        c1 * p.quad_x[j]+c2) 

       p.spline_val .+= p.spline_val_more .* p.quad_w[j] .* c1

   end

   p.spline_val .*= sign * p.delta_x[1]
   
   ind = 1
   for i_grid = index - p.spline_degree:index
      i_mod = mod(i_grid, n_cells ) + 1
      j_dofs[i_mod] += marker_charge * p.spline_val[ind] * p.scaling
      vi = vi - qoverm * p.spline_val[ind] * bfield_dofs[i_mod]
      ind = ind + 1
   end

end


"""
Evaluate field at at position \a position using horner scheme
- `p` : Kernel smoother object 
- `position(p.dim)` : Position of the particle
- `field_dofs_pp(:,:)` : Degrees of freedom in kernel representation.
- `field_value` : Value(s) of the electric fields at given position

""" 
function evaluate_pp(p, position, field_dofs_pp)

    xi = (position[1] - p.domain[1])/p.delta_x[1]
    index = floor(Int64, xi)+1
    xi = xi - (index-1)
   
    horner_1d(p.spline_degree, field_dofs_pp, xi, index)

end

"""
Evaluate field at at position \a position
- `p` : Kernel smoother object 
- `position(p.dim)` : Position of the particle
- `field_dofs(p.n_dofs)` : Coefficient vector for the field DoFs
- `field_value` : Value(s) of the electric fields at given position
"""
function evaluate(p, position, field_dofs)

    xi = (position[1] - p.domain[1])/p.delta_x[1]
    index = floor(Int64, xi)+1
    xi = xi - (index-1)
    index = index - p.spline_degree
    p.spline_val .= uniform_bsplines_eval_basis(p.spline_degree, xi)

    field_value = 0.0
    for i = 1:p.n_span
       index1d = mod(index+i-2, p.n_grid[1])+1
       field_value += field_dofs[index1d] * p.spline_val[i]
    end
    field_value

end

#=

  !---------------------------------------------------------------------------!
  !> Evaluate several fields at position \a position
  subroutine evaluate_multiple_spline_1d(p, position, components, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_1d), intent( inout ) :: p !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(p.dim) !< Position of the particle
    sll_int32,                               intent(in)      :: components(:) !< Components of field_dofs that shall be updated
    sll_real64,                              intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position
    
    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, index
    sll_real64 :: xi[1]

    SLL_ASSERT( size(field_dofs,1) == p.n_dofs )
    SLL_ASSERT( size(field_dofs,2) == size(field_value) )

    xi[1] = (position[1] - p.domain(1))/p.delta_x[1]
    index = ceiling(xi[1])
    xi[1] = xi[1] - (index-1)
    index = index - p.spline_degree    
    call sll_s_uniform_bsplines_eval_basis(p.spline_degree, xi[1], p.spline_val)
    !p.spline_val = sll_f_uniform_b_splines_at_x(p.spline_degree, xi[1])

    field_value = 0.0
    do i1 = 1, p.n_span
       index1d = modulo(index+i1-2, p.n_grid[1])+1
       field_value = field_value + 
            field_dofs(index1d,components) *  
            p.spline_val(i1)
    end do

  end subroutine evaluate_multiple_spline_1d


  

end module sll_m_particle_mesh_coupling_spline_1d
=#

export ParticleMeshCoupling1D

"""
    ParticleMeshCoupling1D( mesh, no_particles, spline_degree, 
                          smoothing_type )
    
Kernel smoother with splines of arbitrary degree placed on a uniform mesh.
Spline with index i starts at point i

- `delta_x` : Value of grid spacing along both directions.
- `n_grid` : Array containing number ofpoints along each direction
- `no_particles` : Number of particles of underlying PIC method 
- `spline_degree` : Degree of smoothing kernel spline
- `n_span` : Number of intervals where spline non zero (spline_degree + 1)
- `scaling` : Scaling factor depending on whether :galerkin or :collocation
- `n_quad_points` : Number of quadrature points
- `spline_val`: scratch data for spline evaluation
- `spline_val_more` : more scratch data for spline evaluation
- `quad_x, quad_w` : quadrature weights and points

!!! note

    Only 1D version is implemented for now

"""
struct ParticleMeshCoupling1D <: AbstractParticleMeshCoupling
    dims::Int
    delta_x::Float64
    n_grid::Int
    n_dofs::Int
    no_particles::Int
    spline_degree::Int
    n_span::Int
    scaling::Float64
    n_quad_points::Int
    spline_val::Vector{Float64}
    spline_val_more::Vector{Float64}
    quad_x::Vector{Float64}
    quad_w::Vector{Float64}
    spline_pp::SplinePP
    xmin::Float64
    Lx::Float64

    function ParticleMeshCoupling1D(
        mesh::OneDGrid, no_particles::Int, spline_degree::Int, smoothing_type::Symbol
    )
        dims = 1
        n_grid = mesh.nx
        n_dofs = n_grid
        xmin = mesh.xmin
        xmax = mesh.xmax
        Lx = xmax - xmin
        delta_x = (xmax - xmin) / n_grid
        n_span = spline_degree + 1

        if smoothing_type == :collocation
            scaling = 1.0 / delta_x
        elseif smoothing_type == :galerkin
            scaling = 1.0
        else
            throw(ArgumentError("""
                 Smoothing Type $smoothing_type not implemented 
                 for kernel_smoother_spline_1d.
           """))
        end

        n_quad_points = (spline_degree + 2) ÷ 2

        spline_val = zeros(Float64, n_span)
        spline_val_more = zeros(Float64, n_span)

        quad_x, quad_w = gausslegendre(n_quad_points)

        spline_pp = SplinePP(spline_degree, n_grid[1])

        return new(
            dims,
            delta_x,
            n_grid,
            n_dofs,
            no_particles,
            spline_degree,
            n_span,
            scaling,
            n_quad_points,
            spline_val,
            spline_val_more,
            quad_x,
            quad_w,
            spline_pp,
            xmin,
            Lx,
        )
    end
end

"""
    add_charge_pp!(rho_dofs, p, position, marker_charge)

Add charge of one particle
- `p`             : kernel smoother object
- `position`      : Position of the particle
- `marker_charge` : Particle weights time charge
- `rho_dofs`      : Coefficient vector of the charge distribution
"""
function add_charge_pp!(
    rho_dofs::Vector{Float64}, p::ParticleMeshCoupling1D, position::Float64, marker_charge
)
    xi = (position - p.xmin) / p.delta_x
    index = trunc(Int, xi) + 1
    xi = xi - (index - 1)
    index = index - p.spline_degree

    @inbounds for i in eachindex(p.spline_val)
        p.spline_val[i] = horner_1d(p.spline_degree, p.spline_pp.poly_coeffs, xi, i)
    end

    @inbounds for i in 1:(p.n_span)
        index1d = mod(index + i - 2, p.n_grid[1]) + 1
        rho_dofs[index1d] += marker_charge * p.spline_val[i] * p.scaling
    end
end

"""  
    add_current_update_v_pp!( j_dofs, p, position_old, position_new, 
                              marker_charge, qoverm, bfield_dofs, vi)

Add current for one particle and update v 
(according to `H_p1` part in Hamiltonian splitting)
"""
function add_current_update_v_pp!(
    j_dofs::AbstractArray,
    p::ParticleMeshCoupling1D,
    position_old,
    position_new,
    marker_charge::Float64,
    qoverm::Float64,
    bfield_dofs::Vector{Float64},
    vi::Vector{Float64},
)

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
        vi[2] = update_jv_pp!(
            j_dofs, p, r_old, r_new, index_old, marker_charge, qoverm, vi[2], bfield_dofs
        )

    elseif index_old < index_new
        vi[2] = update_jv_pp!(
            j_dofs, p, r_old, 1.0, index_old, marker_charge, qoverm, vi[2], bfield_dofs
        )

        vi[2] = update_jv_pp!(
            j_dofs, p, 0.0, r_new, index_new, marker_charge, qoverm, vi[2], bfield_dofs
        )

        for ind in (index_old + 1):(index_new - 1)
            vi[2] = update_jv_pp!(
                j_dofs, p, 0.0, 1.0, ind, marker_charge, qoverm, vi[2], bfield_dofs
            )
        end

    else
        vi[2] = update_jv_pp!(
            j_dofs, p, 1.0, r_new, index_new, marker_charge, qoverm, vi[2], bfield_dofs
        )

        vi[2] = update_jv_pp!(
            j_dofs, p, r_old, 0.0, index_old, marker_charge, qoverm, vi[2], bfield_dofs
        )

        for ind in (index_new + 1):(index_old - 1)
            vi[2] = update_jv_pp!(
                j_dofs, p, 1.0, 0.0, ind, marker_charge, qoverm, vi[2], bfield_dofs
            )
        end
    end
end

"""
    update_jv_pp!( j_dofs, p, lower, upper, index, marker_charge, 
                   qoverm, vi, bfield_dofs)

Helper function for `add_current_update_v`.
"""
function update_jv_pp!(
    j_dofs::AbstractArray,
    p::ParticleMeshCoupling1D,
    lower::Float64,
    upper::Float64,
    index::Int64,
    marker_charge::Float64,
    qoverm::Float64,
    vi::Float64,
    bfield_dofs::Vector{Float64},
)
    n_cells = p.n_grid[1]

    # Evaluation of the primitive integral at the lower 
    # and upper bound of the gridcell

    horner_primitive_1d(p.spline_val, p.spline_degree, p.spline_pp.poly_coeffs_fp, lower)
    horner_primitive_1d(
        p.spline_val_more, p.spline_degree, p.spline_pp.poly_coeffs_fp, upper
    )

    p.spline_val .= (p.spline_val_more .- p.spline_val) .* (p.delta_x[1])

    ind = 1
    @inbounds for i_grid in (index - p.spline_degree):index
        i_mod = mod(i_grid, n_cells) + 1
        j_dofs[i_mod] = j_dofs[i_mod] + (marker_charge * p.spline_val[ind] * p.scaling)
        vi = vi - qoverm * p.spline_val[ind] * bfield_dofs[i_mod]
        ind = ind + 1
    end

    return vi
end

"""
    evaluate_pp(p, position, field_dofs_pp)

Evaluate field at `position` using horner scheme
- `p` : Kernel smoother object 
- `position` : Position of the particle
- `field_dofs_pp` : Degrees of freedom in kernel representation.
- `field_value` : Value(s) of the electric fields at given position

"""
function evaluate_pp(
    p::ParticleMeshCoupling1D, position::Float64, field_dofs_pp::Array{Float64,2}
)
    xi = (position - p.xmin) / p.delta_x
    index = trunc(Int, xi) + 1
    xi = xi - (index - 1)

    return horner_1d(p.spline_degree, field_dofs_pp, xi, index)
end

"""
    add_charge!( rho, p, position, marker_charge) 

Add charge of one particle
- `p`             : kernel smoother object
- `position`      : Position of the particle
- `marker_charge` : Particle weights time charge
- `rho_dofs`      : Coefficient vector of the charge distribution
"""
function add_charge!(
    rho_dofs::Vector{Float64},
    p::ParticleMeshCoupling1D,
    position::Float64,
    marker_charge::Float64,
)
    xi = (position - p.xmin) / p.delta_x
    index = trunc(Int, xi)
    xi = xi - index
    index = index - p.spline_degree

    uniform_bsplines_eval_basis!(p.spline_val, p.spline_degree, xi)

    nx = p.n_grid[1]

    @inbounds for i in 1:(p.n_span)
        index1d = mod1(index + i, nx)
        rho_dofs[index1d] += marker_charge * p.spline_val[i] * p.scaling
    end
end

"""
    add_current_update_v!( j_dofs, p, 
                           position_old, position_new, 
                           marker_charge, qoverm, 
                           bfield_dofs, vi) 

Add current for one particle and update v (according to ``H_{p1}``
part in Hamiltonian splitting)

- Read out particle position and velocity
- Compute index_old, the index of the last DoF on the grid the 
particle contributes to, and `r_old`, its position (normalized to cell size one).

"""
function add_current_update_v!(
    j_dofs::AbstractArray,
    p::ParticleMeshCoupling1D,
    position_old::Float64,
    position_new::Float64,
    marker_charge::Float64,
    qoverm::Float64,
    bfield_dofs::Vector{Float64},
    vi::Float64,
)
    xi = (position_old - p.xmin) / p.delta_x
    index_old = trunc(Int, xi)
    r_old = xi - index_old

    # Compute the new box index index_new and normalized position r_old.

    xi = (position_new - p.xmin) / p.delta_x
    index_new = trunc(Int, xi)
    r_new = xi - index_new

    if index_old == index_new
        if r_old < r_new
            vi = update_jv!(
                j_dofs,
                p,
                r_old,
                r_new,
                index_old,
                marker_charge,
                qoverm,
                1.0,
                vi,
                bfield_dofs,
            )
        else
            vi = update_jv!(
                j_dofs,
                p,
                r_new,
                r_old,
                index_old,
                marker_charge,
                qoverm,
                -1.0,
                vi,
                bfield_dofs,
            )
        end

    elseif index_old < index_new
        vi = update_jv!(
            j_dofs, p, r_old, 1.0, index_old, marker_charge, qoverm, 1.0, vi, bfield_dofs
        )

        vi = update_jv!(
            j_dofs, p, 0.0, r_new, index_new, marker_charge, qoverm, 1.0, vi, bfield_dofs
        )

        for ind in (index_old + 1):(index_new - 1)
            vi = update_jv!(
                j_dofs, p, 0.0, 1.0, ind, marker_charge, qoverm, 1.0, vi, bfield_dofs
            )
        end

    else
        vi = update_jv!(
            j_dofs, p, r_new, 1.0, index_new, marker_charge, qoverm, -1.0, vi, bfield_dofs
        )
        vi = update_jv!(
            j_dofs, p, 0.0, r_old, index_old, marker_charge, qoverm, -1.0, vi, bfield_dofs
        )

        for ind in (index_new + 1):(index_old - 1)
            vi = update_jv!(
                j_dofs, p, 0.0, 1.0, ind, marker_charge, qoverm, -1.0, vi, bfield_dofs
            )
        end
    end

    return vi
end

"""
    update_jv!(j_dofs, p, 
               lower, upper, index, marker_charge, 
               qoverm, sign, vi, bfield_dofs)

Helper function for `add_current_update_v`.
"""
function update_jv!(
    j_dofs::AbstractArray,
    p::ParticleMeshCoupling1D,
    lower::Float64,
    upper::Float64,
    index::Int64,
    marker_charge::Float64,
    qoverm::Float64,
    sign::Float64,
    vi::Float64,
    bfield_dofs::Vector{Float64},
)
    n_cells = p.n_grid[1]

    c1 = 0.5 * (upper - lower)
    c2 = 0.5 * (upper + lower)

    uniform_bsplines_eval_basis!(p.spline_val, p.spline_degree, c1 * p.quad_x[1] + c2)

    p.spline_val .*= p.quad_w[1] * c1

    @inbounds for j in 2:(p.n_quad_points)
        uniform_bsplines_eval_basis!(
            p.spline_val_more, p.spline_degree, c1 * p.quad_x[j] + c2
        )

        p.spline_val .+= p.spline_val_more .* p.quad_w[j] .* c1
    end

    p.spline_val .*= sign * p.delta_x[1]

    ind = 1
    @inbounds for i_grid in (index - p.spline_degree):index
        i_mod = mod(i_grid, n_cells) + 1
        j_dofs[i_mod] += marker_charge * p.spline_val[ind] * p.scaling
        vi = vi - qoverm * p.spline_val[ind] * bfield_dofs[i_mod]
        ind = ind + 1
    end

    return vi
end

"""
    evaluate(p, position, field_dofs)

Evaluate field at `position`

- `p` : Kernel smoother object 
- `position` : Position of the particle
- `field_dofs` : Coefficient vector for the field DoFs
- `field_value` : Value(s) of the electric fields at given position

"""
function evaluate(p::ParticleMeshCoupling1D, position::Float64, field_dofs::Vector{Float64})
    xi = (position - p.xmin) / p.delta_x
    index = trunc(Int, xi)
    dxi = xi - index
    index = index - p.spline_degree
    uniform_bsplines_eval_basis!(p.spline_val, p.spline_degree, dxi)

    nx = p.n_grid[1]
    field_value = 0.0
    @inbounds for i in 1:(p.n_span)
        index1d = mod1(index + i, nx)
        field_value += field_dofs[index1d] * p.spline_val[i]
    end

    return field_value
end

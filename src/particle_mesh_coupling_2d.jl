export ParticleMeshCoupling2D

"""
    ParticleMeshCoupling2D( pg, grid, degree, smoothing_type)

- n_grid(2) : no. of spline coefficients
- domain(2,2) : lower and upper bounds of the domain
- no_particles : no. of particles
- degree : Degree of smoothing kernel spline
- smoothing_type : Define if Galerkin or collocation smoothing for right scaling in accumulation routines 
"""
struct ParticleMeshCoupling2D <: AbstractParticleMeshCoupling

    grid::TwoDGrid
    npart::Int
    spline1::SplinePP
    spline2::SplinePP
    n_span::Int
    degree::Int
    scaling::Float64
    values::Array{Float64,2}

    function ParticleMeshCoupling2D(
        pg::ParticleGroup{2,2},
        grid::TwoDGrid,
        degree::Int,
        smoothing_type::Symbol,
    )

        npart = pg.n_particles
        spline1 = SplinePP(degree, grid.nx)
        spline2 = SplinePP(degree, grid.ny)

        n_span = degree + 1

        if smoothing_type == :collocation
            scaling = 1.0 / (grid.dx * grid.dy)
        elseif smoothing_type == :galerkin
            scaling = 1.0
        else
            println(
                "Smoothing Type $smoothing_type not implemented for kernel_smoother_spline_2d. ",
            )
        end

        values = zeros(n_span, 2)

        new(grid, npart, spline1, spline2, n_span, degree, scaling, values)

    end

end


"""
    compute_shape_factor(pm, xp, yp)

Helper function computing shape factor
- pm : kernel smoother object
"""
function compute_shape_factor(pm::ParticleMeshCoupling2D, xp, yp)

    xp = (xp - pm.grid.xmin) / pm.grid.dx
    yp = (yp - pm.grid.ymin) / pm.grid.dy
    ip = ceil(Int, xp)
    jp = ceil(Int, yp)
    dxp = xp - (ip - 1)
    dyp = yp - (jp - 1)

    uniform_bsplines_eval_basis!(pm.values, pm.degree, dxp, dyp)

    return (ip - pm.degree, jp - pm.degree)

end

"""
    index_1dto2d_column_major(pm, index1d_1, index_1d_2) 

Self function computes the index of a 1D array that stores 2D data in column major ordering. 
It also takes periodic boundary conditions into account.
- index1d_1 !< indice along x (start counting with zero).
- index1d_2 !< indice along y (start counting with zero).
- index2d   !< Corresponding index in 1d array representing 2d data (start counting with one).
"""
function index_1dto2d_column_major(pm, index1d_1, index1d_2)

    index1d_1 = mod(index1d_1, pm.grid.nx)
    index1d_2 = mod(index1d_2, pm.grid.ny)
    index2d = index1d_1 + index1d_2 * pm.grid.nx + 1

    return index2d

end


export add_charge!

"""
    add_charge!(ρ_dofs, pm, xp, yp, wp)

Add charge of single particle
- position : Particle position
- wp : Particle weight times charge
- ρ_dofs : spline coefficient of accumulated density
"""
function add_charge!(ρ_dofs, pm::ParticleMeshCoupling2D, xp, yp, wp)

    ind_x, ind_y = compute_shape_factor(pm, xp, yp)

    @inbounds for i1 = 1:pm.n_span
        index1d_1 = ind_x + i1 - 2
        for i2 = 1:pm.n_span
            index1d_2 = ind_y + i2 - 2
            index2d = index_1dto2d_column_major(pm, index1d_1, index1d_2)
            ρ_dofs[index2d] += (wp * pm.scaling * pm.values[i1, 1] * pm.values[i2, 2])
        end
    end

end

export add_charge_pp!

"""
    add_charge_pp!(ρ_dofs, pm, xp, yp, wp)

## Add charge of single particle

- Information about the 2d mesh
  * delta_x(2)  : Value of grid spacing along both directions.
  *  domain(2,2) : Definition of the domain: domain(1,1) = x1_min, domain(2,1) = x2_min,  domain(1,2) = x1_max, domain(2,2) = x2_max
- Information about the particles
  * no_particles : Number of particles of underlying PIC method (processor local)
  * n_span : Number of intervals where spline non zero (degree + 1)
  * scaling
  
- position : Particle position
- wp : Particle weight times charge
- ρ_dofs : spline coefficient of accumulated density
    
"""
function add_charge_pp!(ρ_dofs, pm::ParticleMeshCoupling2D, xp, yp, wp)

    xp = (xp - pm.grid.xmin) / pm.grid.dx
    yp = (yp - pm.grid.ymin) / pm.grid.dy
    ip = floor(Int, xp) + 1
    jp = floor(Int, yp) + 1
    dxp = xp - (ip - 1)
    dyp = yp - (jp - 1)

    ip = ip - pm.degree
    jp = jp - pm.degree

    horner_m_2d!(pm.values, pm.spline1, pm.spline2, pm.degree, dxp, dyp)

    for i1 = 1:pm.n_span
        index1d_1 = ip + i1 - 2
        for i2 = 1:pm.n_span
            index1d_2 = jp + i2 - 2
            index2d = index_1dto2d_column_major(pm, index1d_1, index1d_2)
            ρ_dofs[index2d] += (wp * pm.scaling * pm.values[i1, 1] * pm.values[i2, 2])
        end
    end

end

export evaluate_pp, evaluate

"""
    evaluate_pp(pm, xp, yp, pp)

Evaluate field at position using horner scheme
- pm : kernel smoother object    
- position : Position where to evaluate
- field_dofs_pp(:,:) : Degrees of freedom in kernel representation.
- field_value : Value of the field

"""
function evaluate_pp(pm::ParticleMeshCoupling2D, xp, yp, pp)

    xi = (xp - pm.grid.xmin) / pm.grid.dx
    yi = (yp - pm.grid.ymin) / pm.grid.dy

    idx = floor(Int, xi) + 1
    idy = floor(Int, yi) + 1

    xi = xi - (idx - 1)
    yi = yi - (idy - 1)

    d = pm.degree

    nx = pm.grid.nx
    ny = pm.grid.ny

    horner_2d((d, d), pp, (xi, yi), (idx, idy), (nx, ny))

end


"""
    evaluate(pm, xp, yp, field_dofs)

- position(pm.dim) : Position where to evaluate
- field_dofs(pm.n_dofs) : Degrees of freedom in kernel representation.
- field_value : Value of the field

Evaluate field with given dofs at position
"""
function evaluate(pm, xp, yp, field_dofs)

    ix, iy = compute_shape_factor(pm, xp, yp)

    value = 0.0
    for i1 = 1:pm.n_span
        index1d_1 = ix + i1 - 2
        for i2 = 1:pm.n_span
            index1d_2 = iy + i2 - 2
            index2d = index_1dto2d_column_major(pm, index1d_1, index1d_2)
            value += field_dofs[index2d] * pm.values[i1, 1] * pm.values[i2, 2]
        end
    end

    value

end

"""
    evaluate_multiple(pm, position, field_dofs)

## Evaluate multiple fields at position \a position
- position(pm%dim) : Position where to evaluate
- components(:) : Components of the field that shall be evaluated
- field_dofs(:,:) : Degrees of freedom in kernel representation.
- field_value(:) : Value of the field
"""
function evaluate_multiple(pm, position, field_dofs)

    indices = compute_shape_factor(pm, position...)

    field_value1 = 0.0
    field_value2 = 0.0
    for i1 = 1:pm.n_span
        index1d_1 = indices[1] + i1 - 2
        for i2 = 1:pm.n_span
            index1d_2 = indices[2] + i2 - 2
            index2d = index_1dto2d_column_major(pm, index1d_1, index1d_2)
            c = pm.values[i1, 1] * pm.values[i2, 2]
            field_value1 += field_dofs[1][index2d] * c
            field_value2 += field_dofs[2][index2d] * c
        end
    end

    field_value1, field_value2

end

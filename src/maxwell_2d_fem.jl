export TwoDMaxwell

"""
    TwoDMaxwell( mesh, degree)

- mesh       : cartesian mesh
- s_deg_0    : spline degree 0-forms
- s_deg_1    : spline degree 1-forms

"""
struct TwoDMaxwell
    s_deg_0::Int
    s_deg_1::Int

    mesh::TwoDGrid

    mass_line_0::Vector{Array{Float64,1}}
    mass_line_1::Vector{Array{Float64,1}}
    mass_line_mixed::Vector{Array{Float64,1}}
    inv_mass_1::Vector{TwoDLinearSolverSplineMass}
    inv_mass_2::Vector{TwoDLinearSolverSplineMass}
    poisson::TwoDPoisson
    wk1d::Array{Float64,1}
    wk2d::Array{Float64,2}

    function TwoDMaxwell(mesh, degree)
        nx, ny = mesh.nx, mesh.ny
        dx, dy = mesh.dx, mesh.dy
        s_deg_0 = degree
        s_deg_1 = degree - 1

        # Sparse matrices
        # Assemble the mass matrices
        # First assemble a mass line for both degrees

        mass_line_0 = [
            spline_fem_mass_line(s_deg_0) .* dx, spline_fem_mass_line(s_deg_0) .* dy
        ]

        mass_line_1 = [
            spline_fem_mass_line(s_deg_1) .* dx, spline_fem_mass_line(s_deg_1) .* dy
        ]

        mass_line_mixed = [
            spline_fem_mixedmass_line(s_deg_0) .* dx,
            spline_fem_mixedmass_line(s_deg_0) .* dy,
        ]

        # Next put together the 1d parts of the 2d Kronecker product

        eig_values_mass_0_1 = spline_fem_compute_mass_eig(nx, s_deg_0, mass_line_0[1])
        eig_values_mass_0_2 = spline_fem_compute_mass_eig(ny, s_deg_0, mass_line_0[2])
        eig_values_mass_1_1 = spline_fem_compute_mass_eig(nx, s_deg_1, mass_line_1[1])
        eig_values_mass_1_2 = spline_fem_compute_mass_eig(ny, s_deg_1, mass_line_1[2])

        inv_mass_1 = [
            TwoDLinearSolverSplineMass(nx, ny, eig_values_mass_1_1, eig_values_mass_0_2),
            TwoDLinearSolverSplineMass(nx, ny, eig_values_mass_0_1, eig_values_mass_1_2),
            TwoDLinearSolverSplineMass(nx, ny, eig_values_mass_0_1, eig_values_mass_0_2),
        ]

        inv_mass_2 = [
            TwoDLinearSolverSplineMass(nx, ny, eig_values_mass_0_1, eig_values_mass_1_2),
            TwoDLinearSolverSplineMass(nx, ny, eig_values_mass_1_1, eig_values_mass_0_2),
            TwoDLinearSolverSplineMass(nx, ny, eig_values_mass_1_1, eig_values_mass_1_2),
        ]

        poisson = TwoDPoisson(mesh, s_deg_0)

        wk1d = zeros(nx * ny)
        wk2d = zeros(nx, ny)

        return new(
            s_deg_0,
            s_deg_1,
            mesh,
            mass_line_0,
            mass_line_1,
            mass_line_mixed,
            inv_mass_1,
            inv_mass_2,
            poisson,
            wk1d,
            wk2d,
        )
    end
end

function spline_fem_mixedmass_line(deg)
    n = min(3 * deg + 1, 10)

    spline_val_0 = zeros(deg + 1, n)
    spline_val_1 = zeros(deg, n)

    x, w = gausslegendre(n)
    x .+= 1
    x ./= 2
    w ./= 2

    for j in 1:n
        spline_val_0[:, j] .= uniform_bsplines_eval_basis(deg, x[j])
        spline_val_1[:, j] .= uniform_bsplines_eval_basis(deg - 1, x[j])
    end

    mass_line = zeros(2deg)
    for j in 2:(deg + 1), i in j:(deg + 1), k in 1:n
        mass_line[j + deg - 1] += spline_val_0[i, k] * spline_val_1[i - j + 1, k] * w[k]
    end

    for j in (-deg + 1):0, i in 1:(deg + j), k in 1:n
        mass_line[j + deg] += spline_val_0[i, k] * spline_val_1[i - j, k] * w[k]
    end

    return mass_line
end

"""
    compute_rhs_from_function(solver, func, component, form)

Compute the FEM right-hand-side for a given function f and periodic splines of given degree
Its components are ``\\int f N_i dx`` where ``N_i`` is the B-spline starting at ``x_i`` 
"""
function compute_rhs_from_function(
    solver::TwoDMaxwell, f::Function, component::Int, form::Int
)
    nx, ny = solver.mesh.nx, solver.mesh.ny
    dx, dy = solver.mesh.dx, solver.mesh.dy
    deg_0 = solver.s_deg_0
    deg_1 = solver.s_deg_1

    coefs_dofs = zeros(nx * ny)  # Finite Element right-hand-side

    # Define the spline degree in the 3 dimensions, depending on form and component of the form
    if form == 0
        degree = [deg_0, deg_0]
    elseif (form == 1)
        degree = [deg_0, deg_0]
        component < 3 && (degree[component] = deg_1)
    elseif form == 2
        degree = [deg_1, deg_1]
        component < 3 && (degree[component] = deg_0)
    elseif form == 3
        degree = [deg_1, deg_1]
    else
        @error " Wrong form "
    end

    d1, d2 = degree

    # take enough Gauss points so that projection is exact for splines of degree deg
    # rescale on [0,1] for compatibility with B-splines
    bspl_d1 = zeros(d1 + 1, d1 + 1)
    x1, w1 = gausslegendre(d1 + 1)
    x1 .+= 1
    x1 ./= 2
    w1 ./= 2
    # Compute bsplines at gauss_points
    for k in 1:(d1 + 1)
        bspl_d1[k, :] .= uniform_bsplines_eval_basis(d1, x1[k])
    end

    bspl_d2 = zeros(d2 + 1, d2 + 1)
    x2, w2 = gausslegendre(d2 + 1)
    x2 .+= 1
    x2 ./= 2
    w2 ./= 2
    # Compute bsplines at gauss_points
    for k in 1:(d2 + 1)
        bspl_d2[k, :] .= uniform_bsplines_eval_basis(d2, x2[k])
    end

    counter = 1
    # Compute coefs_dofs = int f(x)N_i(x) 
    for i2 in 1:ny, i1 in 1:nx
        coef = 0.0
        # loop over support of B spline
        for j1 in 1:(d1 + 1), j2 in 1:(d2 + 1)
            # loop over Gauss points
            for k1 in 1:(d1 + 1), k2 in 1:(d2 + 1)
                x = dx * (x1[k1] + i1 + j1 - 2)
                y = dy * (x2[k2] + i2 + j2 - 2)
                coef +=
                    w1[k1] *
                    w2[k2] *
                    f(x, y) *
                    bspl_d1[k1, d1 + 2 - j1] *
                    bspl_d2[k2, d2 + 2 - j2]
            end
        end

        # rescale by cell size
        coefs_dofs[counter] = coef * dx * dy
        counter = counter + 1
    end

    return coefs_dofs
end

function compute_e_from_rho!(efield, solver::TwoDMaxwell, rho)
    return compute_e_from_rho!(efield, solver.poisson, rho)
end

"""
    compute_fem_rhs(solver, func , component, form, coefs_dofs)

Compute the FEM right-hand-side for a given function f and periodic splines of given degree
Its components are ``\\int f N_i dx`` where ``N_i`` is the B-spline starting at ``x_i`` 
- component : Specify the component
- form : Specify 0,1,2 or 3 form
coefs_dofs - Finite Element right-hand-side
"""
function compute_fem_rhs!(
    coefs_dofs::Vector{Float64}, solver::TwoDMaxwell, f, component::Int, form::Int
)
    n1, n2 = solver.mesh.nx, solver.mesh.ny
    δ1, δ2 = solver.mesh.dx, solver.mesh.dy
    s_deg_0, s_deg_1 = solver.s_deg_0, solver.s_deg_1

    if form == 0
        degree = [s_deg_0, s_deg_0]
    elseif form == 1
        degree = [s_deg_0, s_deg_0]
        if component < 3
            degree[component] = s_deg_1
        end
    elseif form == 2
        degree = [s_deg_1, s_deg_1]
        if component < 3
            degree[component] = s_deg_0
        end
    elseif form == 3
        degree = [s_deg_1, s_deg_1]
    end

    d1, d2 = degree

    # take enough Gauss points so that projection is exact for splines of degree deg
    # rescale on [0,1] for compatibility with B-splines

    bspl_d1 = zeros(d1 + 1, d1 + 1)
    x1, w1 = gausslegendre(d1 + 1)
    x1 .+= 1
    x1 ./= 2
    w1 ./= 2
    # Compute bsplines at gauss_points
    for k in 1:(d1 + 1)
        bspl_d1[k, :] .= uniform_bsplines_eval_basis(d1, x1[k])
    end

    bspl_d2 = zeros(d2 + 1, d2 + 1)
    x2, w2 = gausslegendre(d2 + 1)
    x2 .+= 1
    x2 ./= 2
    w2 ./= 2
    # Compute bsplines at gauss_points
    for k in 1:(d2 + 1)
        bspl_d2[k, :] .= uniform_bsplines_eval_basis(d2, x2[k])
    end

    # Compute coefs_dofs = int f(x)N_i(x) 
    ind = 0
    for i2 in 1:n2, i1 in 1:n1
        ind += 1
        coef = 0.0
        # loop over support of B spline
        for j1 in 1:(d1 + 1), j2 in 1:(d2 + 1)
            # loop over Gauss points
            for k1 in 1:(d1 + 1), k2 in 1:(d2 + 1)
                xg1 = δ1 * (x1[k1] + i1 + j1 - 2)
                xg2 = δ2 * (x2[k2] + i2 + j2 - 2)
                coef +=
                    w1[k1] *
                    w2[k2] *
                    f(xg1, xg2) *
                    bspl_d1[k1, d1 + 2 - j1] *
                    bspl_d2[k2, d2 + 2 - j2]
            end
        end
        # rescale by cell size
        coefs_dofs[ind] = coef * δ1 * δ2
    end
end

export l2projection
"""
Compute the L2 projection of a given function f on periodic splines of given degree
"""
function l2projection(solver::TwoDMaxwell, f, component, form)
    n1, n2 = solver.mesh.nx, solver.mesh.ny

    compute_fem_rhs!(solver.wk1d, solver, f, component, form)

    if form == 1
        solve(solver.inv_mass_1[component], solver.wk1d)
    elseif form == 2
        solve(solver.inv_mass_2[component], solver.wk1d)
    end
end

function spline_fem_multiply_mass(n_cells, degree, mass, invec)
    outvec = copy(invec)

    ind = 1
    # For the first degree rows we need to put the first part to the back due to periodic boundaries
    for row in 1:degree
        outvec[row] = mass[1] * invec[row]
        for column in 1:(row - 1)
            outvec[row] += mass[column + 1] * (invec[row + column] + invec[row - column])
        end
        for column in row:degree
            outvec[row] +=
                mass[column + 1] * (invec[row + column] + invec[row - column + n_cells])
        end
    end

    for row in (degree + 1):(n_cells - degree)
        outvec[row] = mass[1] * invec[row]
        for column in 1:degree
            outvec[row] += mass[column + 1] * (invec[row + column] + invec[row - column])
        end
    end

    # For the last degree rows, we need to put the second part to the front due to periodic boundaries
    for row in (n_cells - degree + 1):n_cells
        outvec[row] = mass[1] * invec[row]
        for column in 1:(n_cells - row)
            outvec[row] += mass[column + 1] * (invec[row + column] + invec[row - column])
        end
        for column in (n_cells - row + 1):degree
            outvec[row] +=
                mass[column + 1] * (invec[row + column - n_cells] + invec[row - column])
        end
    end

    return outvec
end

"""
    multiply_mass_2dkron!( c_out, solver, mass_line_1, mass_line_2,  c_in )

Multiply by the mass matrix 
"""
function multiply_mass_2dkron!(c_out, solver, mass_line_1, mass_line_2, c_in)
    nx1 = solver.mesh.nx
    nx2 = solver.mesh.ny
    deg1 = size(mass_line_1)[1] - 1
    deg2 = size(mass_line_2)[1] - 1

    solver.wk2d .= reshape(c_in, nx1, nx2)

    for j in 1:nx2
        solver.wk2d[:, j] .= spline_fem_multiply_mass(
            nx1, deg1, mass_line_1, solver.wk2d[:, j]
        )
    end
    for i in 1:nx1
        solver.wk2d[i, :] .= spline_fem_multiply_mass(
            nx2, deg2, mass_line_2, solver.wk2d[i, :]
        )
    end

    return c_out .= vec(solver.wk2d)
end

"""
    compute_e_from_b!(e, solver, delta_t, b)

compute Ey from Bz using weak Ampere formulation 
"""
function compute_e_from_b!(e, solver, dt, b)
    nx1, nx2 = solver.mesh.nx, solver.mesh.ny
    dx1, dx2 = solver.mesh.dx, solver.mesh.dy

    work = deepcopy(b)

    multiply_mass_2dkron!(
        work[1], solver, solver.mass_line_0[1], solver.mass_line_1[2], b[1]
    )
    multiply_mass_2dkron!(
        work[2], solver, solver.mass_line_1[1], solver.mass_line_0[2], b[2]
    )
    multiply_mass_2dkron!(
        work[3], solver, solver.mass_line_1[1], solver.mass_line_1[2], b[3]
    )

    curl_b = deepcopy(work)

    # Compute curl with periodic boundary conditions
    for j in 1:nx2
        indp2 = j == nx2 ? nx1 * (1 - nx2) : nx1
        for i in 1:nx1
            indp1 = i == nx1 ? 1 - nx1 : 1
            ind2d = (j - 1) * nx1 + i
            ind2d_1 = ind2d + indp1
            ind2d_2 = ind2d + indp2
            curl_b[1][ind2d] = -(work[3][ind2d] - work[3][ind2d_2]) / dx2
            curl_b[2][ind2d] = (work[3][ind2d] - work[3][ind2d_1]) / dx1
            curl_b[3][ind2d] =
                (work[1][ind2d] - work[1][ind2d_2]) / dx2 -
                (work[2][ind2d] - work[2][ind2d_1]) / dx1
        end
    end

    curl_b[1] .= solve(solver.inv_mass_1[1], curl_b[1])
    curl_b[2] .= solve(solver.inv_mass_1[2], curl_b[2])
    curl_b[3] .= solve(solver.inv_mass_1[3], curl_b[3])

    # Update b from solver value
    e[1] .= e[1] .+ dt .* curl_b[1]
    e[2] .= e[2] .+ dt .* curl_b[2]
    return e[3] .= e[3] .+ dt .* curl_b[3]
end

"""
    compute_b_from_e!(b, solver, delta_t, e)

Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
``B_z^{new}(x_j) = B_z^{old}(x_j) - \\frac{\\Delta t}{\\Delta x} (E_y(x_j) - E_y(x_{j-1}) ``

multiplying by discrete curl matrix

"""
function compute_b_from_e!(b, solver, dt, e)
    nx1 = solver.mesh.nx
    nx2 = solver.mesh.ny
    dx1 = solver.mesh.dx
    dx2 = solver.mesh.dy

    for j in 1:nx2
        indp2 = j == 1 ? nx1 * (nx2 - 1) : -nx1
        for i in 1:nx1
            indp1 = i == 1 ? nx1 - 1 : -1
            ind2d = (j - 1) * nx1 + i
            ind2d_1 = ind2d + indp1
            ind2d_2 = ind2d + indp2

            b[1][ind2d] += -dt * (e[3][ind2d] - e[3][ind2d_2]) / dx2
            b[2][ind2d] += dt * (e[3][ind2d] - e[3][ind2d_1]) / dx1
            b[3][ind2d] +=
                -dt *
                ((e[2][ind2d] - e[2][ind2d_1]) / dx1 - (e[1][ind2d] - e[1][ind2d_2]) / dx2)
        end
    end
end

"""
    compute_e_from_j(e, solver, current, component)

Compute E_i from j_i integrated over the time interval using weak Ampere formulation
- solver : Maxwell solver class
- current : Component of the current integrated over time interval
- component : Component of the Efield to be computed
- e : Updated electric field
"""
function compute_e_from_j!(e, solver::TwoDMaxwell, current, component)
    work = solve(solver.inv_mass_1[component], current)

    return e .= e .- work
end

export compute_rho_from_e!

"""
    compute_rho_from_e(rho, solver, efield)

compute rho from e using weak Gauss law ( rho = G^T M_1 e ) 
"""
function compute_rho_from_e!(rho, solver, efield)
    work = deepcopy(efield)
    multiply_mass_2dkron!(
        work[1], solver, solver.mass_line_1[1], solver.mass_line_0[2], efield[1]
    )
    multiply_mass_2dkron!(
        work[2], solver, solver.mass_line_0[1], solver.mass_line_1[2], efield[2]
    )
    multiply_mass_2dkron!(
        work[3], solver, solver.mass_line_0[1], solver.mass_line_0[2], efield[3]
    )

    dx1, dx2 = solver.mesh.dx, solver.mesh.dy
    nx1, nx2 = solver.mesh.nx, solver.mesh.ny

    # Compute div with periodic boundary conditions
    for j in 1:nx2
        indp2 = j == nx2 ? -nx1 * (nx2 - 1) : nx1
        for i in 1:nx1
            indp1 = i == nx1 ? 1 - nx1 : 1
            ind2d = (j - 1) * nx1 + i
            ind2d_1 = ind2d + indp1
            ind2d_2 = ind2d + indp2
            rho[ind2d] = (
                (work[1][ind2d] - work[1][ind2d_1]) / dx1 +
                (work[2][ind2d] - work[2][ind2d_2]) / dx2
            )
        end
    end

    return rho .*= -1
end

"""
    inner_product( solver, coefs1_dofs, coefs2_dofs, component, form)

- solver : Maxwell solver object
- coefs1_dofs : Coefficient for each DoF
- coefs2_dofs : Coefficient for each DoF
- component : Specify the component
- form : Specify 0,1,2 or 3 form
- r : Result: squared L2 norm
"""
function inner_product(solver, coefs1_dofs, coefs2_dofs, component, form)
    if form == 0
        multiply_mass_2dkron!(
            solver.wk1d, solver, solver.mass_line_0[1], solver.mass_line_0[2], coefs2_dofs
        )
    elseif form == 1
        if component == 1
            multiply_mass_2dkron!(
                solver.wk1d,
                solver,
                solver.mass_line_1[1],
                solver.mass_line_0[2],
                coefs2_dofs,
            )
        elseif component == 2
            multiply_mass_2dkron!(
                solver.wk1d,
                solver,
                solver.mass_line_0[1],
                solver.mass_line_1[2],
                coefs2_dofs,
            )
        elseif component == 3
            multiply_mass_2dkron!(
                solver.wk1d,
                solver,
                solver.mass_line_0[1],
                solver.mass_line_0[2],
                coefs2_dofs,
            )
        end
    elseif form == 2
        if component == 1
            multiply_mass_2dkron!(
                solver.wk1d,
                solver,
                solver.mass_line_0[1],
                solver.mass_line_1[2],
                coefs2_dofs,
            )
        elseif component == 2
            multiply_mass_2dkron!(
                solver.wk1d,
                solver,
                solver.mass_line_1[1],
                solver.mass_line_0[2],
                coefs2_dofs,
            )
        elseif component == 3
            multiply_mass_2dkron!(
                solver.wk1d,
                solver,
                solver.mass_line_1[1],
                solver.mass_line_1[2],
                coefs2_dofs,
            )
        end
    elseif form == 3
        multiply_mass_2dkron!(
            solver.wk1d, solver, solver.mass_line_1[1], solver.mass_line_1[2], coefs2_dofs
        )
    end

    return sum(coefs1_dofs .* solver.wk1d)
end

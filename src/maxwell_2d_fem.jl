export TwoDMaxwell

"""
    TwoDMaxwell( mesh, degree)

- mesh       : cartesian mesh
- s_deg_0    : spline degree 0-forms
- s_deg_1    : spline degree 1-forms

"""
struct TwoDMaxwell

    s_deg_0 :: Int
    s_deg_1 :: Int

    mesh :: TwoDGrid

    mass_line_0_x :: Array{Float64,1}
    mass_line_1_x :: Array{Float64,1}
    mass_line_mixed_x :: Array{Float64,1}
    mass_line_0_y :: Array{Float64,1}
    mass_line_1_y :: Array{Float64,1}
    mass_line_mixed_y :: Array{Float64,1}

    inverse_mass_1 :: Vector{TwoDLinearSolverSplineMass}
    inverse_mass_2 :: Vector{TwoDLinearSolverSplineMass}
    poisson :: TwoDPoisson

    function TwoDMaxwell( mesh, degree )

        nx, ny = mesh.nx, mesh.ny
        dx, dy = mesh.dx, mesh.dy
        s_deg_0 = degree
        s_deg_1 = degree - 1

        # Sparse matrices
        # Assemble the mass matrices
        # First assemble a mass line for both degrees

        mass_line_0_x = spline_fem_mass_line( s_deg_0 ) .* dx
        mass_line_0_y = spline_fem_mass_line( s_deg_0 ) .* dy
        mass_line_1_x = spline_fem_mass_line( s_deg_1 ) .* dx
        mass_line_1_y = spline_fem_mass_line( s_deg_1 ) .* dy
        mass_line_mixed_x = spline_fem_mixedmass_line( s_deg_0 ) .* dx
        mass_line_mixed_y = spline_fem_mixedmass_line( s_deg_0 ) .* dy

        # Next put together the 1d parts of the 2d Kronecker product

        eig_values_mass_0_1 = spline_fem_compute_mass_eig( nx, s_deg_0, mass_line_0_x)
        eig_values_mass_0_2 = spline_fem_compute_mass_eig( ny, s_deg_0, mass_line_0_y)
        eig_values_mass_1_1 = spline_fem_compute_mass_eig( nx, s_deg_1, mass_line_1_x)
        eig_values_mass_1_2 = spline_fem_compute_mass_eig( ny, s_deg_1, mass_line_1_y)

        inverse_mass_1 = TwoDLinearSolverSplineMass[]
        push!(inverse_mass_1, TwoDLinearSolverSplineMass( nx, ny, eig_values_mass_1_1, eig_values_mass_0_2))
        push!(inverse_mass_1, TwoDLinearSolverSplineMass( nx, ny, eig_values_mass_0_1, eig_values_mass_1_2))
        push!(inverse_mass_1, TwoDLinearSolverSplineMass( nx, ny, eig_values_mass_0_1, eig_values_mass_0_2))

        inverse_mass_2 = TwoDLinearSolverSplineMass[]
        push!(inverse_mass_2, TwoDLinearSolverSplineMass( nx, ny, eig_values_mass_0_1, eig_values_mass_1_2))
        push!(inverse_mass_2, TwoDLinearSolverSplineMass( nx, ny, eig_values_mass_1_1, eig_values_mass_0_2))
        push!(inverse_mass_2, TwoDLinearSolverSplineMass( nx, ny, eig_values_mass_1_1, eig_values_mass_1_2))

        poisson = TwoDPoisson( mesh, s_deg_0 )

        new( s_deg_0, s_deg_1, mesh, 
             mass_line_0_x, mass_line_1_x, mass_line_mixed_x, 
             mass_line_0_y, mass_line_1_y, mass_line_mixed_y,
             inverse_mass_1, inverse_mass_2, poisson )
     
    end 

end 

function spline_fem_mixedmass_line( deg )

    n = min(3*deg+1,10)

    spline_val_0 = zeros(deg+1, n)
    spline_val_1 = zeros(deg, n)

    x, w = gausslegendre( n )
    x .+= 1; x ./= 2; w ./= 2

    for j=1:n
        spline_val_0[:,j] .= uniform_bsplines_eval_basis( deg, x[j])
        spline_val_1[:,j] .= uniform_bsplines_eval_basis( deg-1, x[j])
    end

    mass_line = zeros(2deg)
    for j=2:deg+1, i= j:deg+1, k = 1:n
        mass_line[j+deg-1] += spline_val_0[i,k] * spline_val_1[i-j+1,k]*w[k]
    end

    for j=-deg+1:0, i= 1:deg+j, k = 1:n
        mass_line[j+deg] += spline_val_0[i,k] * spline_val_1[i-j,k]*w[k]
    end

    mass_line


end 


"""
    compute_rhs_from_function(solver, func, component, form)

Compute the FEM right-hand-side for a given function f and periodic splines of given degree
Its components are ``\\int f N_i dx`` where ``N_i`` is the B-spline starting at ``x_i`` 
"""
function compute_rhs_from_function(solver :: TwoDMaxwell, f :: Function, component, form)

     nx, ny = solver.mesh.nx, solver.mesh.ny
     dx, dy = solver.mesh.dx, solver.mesh.dy
     deg_0 = solver.s_deg_0
     deg_1 = solver.s_deg_1

     coefs_dofs = zeros(nx * ny)  # Finite Element right-hand-side

     degree = zeros(Int, 2)
     # Define the spline degree in the 3 dimensions, depending on form and component of the form
     if form == 0
        degree .= deg_0
     elseif (form == 1 )
        degree .= deg_0
        if component<3
           degree[component] = deg_1
        end
     elseif form == 2
        degree .= deg_1
        if (component<3)
           degree[component] = deg_0
        end
     elseif form == 3
        degree .=  deg_1
     else 
        @error " Wrong form "
     end

     d1, d2 = degree

     # take enough Gauss points so that projection is exact for splines of degree deg
     # rescale on [0,1] for compatibility with B-splines
     bspl_d1 = zeros(d1+1, d1+1)
     x1, w1  = gausslegendre(d1+1)
     x1 .+= 1; x1 ./= 2; w1 ./= 2
     # Compute bsplines at gauss_points
     for k=1:d1+1
         bspl_d1[k,:] .= uniform_bsplines_eval_basis(d1,x1[k])
     end

     bspl_d2 = zeros(d2+1, d2+1)
     x2, w2 = gausslegendre(d2+1)
     x2 .+= 1; x2 ./= 2; w2 ./= 2
     # Compute bsplines at gauss_points
     for k=1:d2+1
         bspl_d2[k,:] .= uniform_bsplines_eval_basis(d2,x2[k])
     end

     counter = 1
     # Compute coefs_dofs = int f(x)N_i(x) 
     for i2 = 1:ny, i1 = 1:nx
         coef = 0.0
         # loop over support of B spline
         for j1 = 1:d1+1, j2 = 1:d2+1
             # loop over Gauss points
             for k1=1:d1+1, k2=1:d2+1
                 x = dx*(x1[k1] + i1 + j1 - 2)
                 y = dy*(x2[k2] + i2 + j2 - 2)
                 coef += w1[k1] * w2[k2] * f(x, y) * bspl_d1[k1,d1+2-j1]* bspl_d2[k2,d2+2-j2]
             end
         end

         # rescale by cell size
         coefs_dofs[counter] = coef * dx * dy
         counter = counter+1
     end

end 

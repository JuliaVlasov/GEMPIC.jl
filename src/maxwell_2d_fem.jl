export TwoDMaxwell


struct TwoDMaxwell

    nx :: Int
    ny :: Int
    degree :: Int
    dx :: Float64
    dy :: Float64

    eig_values_dtm1d_1 :: Vector{Float64}  
    eig_values_dtm1d_2 :: Vector{Float64}
    eig_values_d1 :: Vector{ComplexF64}    
    eig_values_d2 :: Vector{ComplexF64}
    eig_values_mass_0_1 :: Vector{Float64}
    eig_values_mass_0_2 :: Vector{Float64}

    array1d_x :: Vector{ComplexF64}
    array1d_y :: Vector{ComplexF64}
    scratch   :: Array{ComplexF64, 2}
    scratchx  :: Array{ComplexF64, 2}
    scratchy  :: Array{ComplexF64, 2}

    function TwoDMaxwell( mesh, degree )

        nx, ny = (mesh.nx, mesh.ny)
        dx, dy = mesh.dx, mesh.dy

        mass_line_0 = zeros(degree+1)
        mass_line_1 = zeros(degree)
        eig_values_mass_1_1 = zeros(nx)
        eig_values_mass_1_2 = zeros(ny)

        array1d_x = zeros(ComplexF64, nx)
        array1d_y = zeros(ComplexF64, ny)
        scratch  = zeros(ComplexF64, nx, ny)
        scratchx = zeros(ComplexF64, nx, ny)
        scratchy = zeros(ComplexF64, nx, ny)

        # Eigenvalues of mass matrices
        mass_line_0 = spline_fem_mass_line( degree )
        mass_line_1 = spline_fem_mass_line( degree-1)
        
        eig_values_mass_0_1 = spline_fem_compute_mass_eig( nx, degree, mass_line_0*dx )
        eig_values_mass_1_1 = spline_fem_compute_mass_eig( nx, degree-1, mass_line_1*dx )

        eig_values_mass_0_2 = spline_fem_compute_mass_eig( ny, degree, mass_line_0*dy )
        eig_values_mass_1_2 = spline_fem_compute_mass_eig( ny, degree-1, mass_line_1*dy )

        eig_values_d1 = zeros(ComplexF64, nx)
        eig_values_dtm1d_1 = zeros(nx)

        for j=2:nx
            angle = 2π*(j-1)/nx
            eig_values_d1[j] = (1 - cos(angle))/dx + 1im * sin(angle)/dx
            eig_values_dtm1d_1[j] = 2/dx^2 * (1-cos(angle))* eig_values_mass_1_1[j]
        end

        eig_values_d2 = zeros(ComplexF64, ny)
        eig_values_dtm1d_2 = zeros(ny)

        for j=2:ny
            angle = 2π*(j-1)/ny
            eig_values_d2[j] = (1 - cos(angle))/dy + 1im * sin(angle)/dy
            eig_values_dtm1d_2[j] = 2/dy^2*(1-cos(angle))* eig_values_mass_1_2[j]
        end

        new( nx, ny, degree, dx, dy, eig_values_dtm1d_1, eig_values_dtm1d_2,
             eig_values_d1, eig_values_d2, 
             eig_values_mass_0_1, eig_values_mass_0_2,
             array1d_x, array1d_y,
             scratch, scratchx, scratchy )

    end
     
end

export compute_rhs_from_function

"""
    compute_rhs_from_function(solver, func, coefs_dofs)

Compute the FEM right-hand-side for a given function f and periodic splines of given degree
Its components are ``\int f N_i dx`` where ``N_i`` is the B-spline starting at ``x_i`` 
- solver :: Maxwell solver 
- func : function
- coefs_dofs[:]  : Finite Element right-hand-side
"""
function compute_rhs_from_function( solver :: TwoDMaxwell, f)

     nx, ny = solver.nx, solver.ny
     dx, dy = solver.dx, solver.dy
     d = solver.degree

     rhs = zeros(nx * ny)

     # take enough Gauss points so that projection is exact for splines of degree deg
     # rescale on [0,1] for compatibility with B-splines
     bspl_d1 = zeros(d+1, d+1)
     xg, wg = gausslegendre(d+1)
     xg .+= 1
     xg ./= 2
     wg ./= 2
     # Compute bsplines at gauss_points
     for k=1:d+1
         bspl_d1[k,:] .= uniform_bsplines_eval_basis(d, xg[k]) 
     end

     counter = 1
     # Compute coefs_dofs = int f(x)N_i(x) 
     for i2 = 1:ny, i1 = 1:nx
         coef=0.0
         # loop over support of B spline
         for j1 = 1:d+1, j2 = 1:d+1
             # loop over Gauss points
             for k1=1:d+1, k2=1:d+1
                 x = dx*(xg[k1] + (i1 + j1 - 2))
                 y = dy*(xg[k2] + (i2 + j2 - 2))
                 coef += wg[k1]* wg[k2] * f(x, y) * bspl_d1[k1,d+2-j1] * bspl_d1[k2,d+2-j2]
             end
         end
         # rescale by cell size
         rhs[counter] = coef * dx * dy
         counter += 1
     end

     rhs

end

"""
Compute the FEM right-hand-side for a given function f and periodic splines of given degree
Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
"""
function compute_rhs_from_function(solver, func, component, form)

     nx, ny = solver.nx, solver.ny
     dx, dy = solver.dx, solver.dy
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
     x1 .+= 1
     x1 ./= 2
     w1 ./= 2
     # Compute bsplines at gauss_points
     for k=1:d1+1
         bspl_d1[k,:] .= uniform_bsplines_eval_basis(d1,x1[k])
     end

     bspl_d2 = zeros(d2+1, d2+1)
     x2, w2 = gausslegendre(d2+1)
     x2 .+= 1
     x2 ./= 2
     w2 ./= 2
     # Compute bsplines at gauss_points
     for k=1,d2+1
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


function spline_fem_mass_line ( degree )


    q = degree+1

    quad_xw = zeros(2, q)
    spline_val = zeros( degree+1, q )

    quad_x, quad_w = gausslegendre( q )
    quad_x .+= 1
    quad_x .*= 0.5
    quad_w .*= 0.5

    for j=1:q
        spline_val[:,j] .= uniform_bsplines_eval_basis( degree, quad_x[j] )
    end

    mass_line = zeros(degree+1)
    for j=1:degree+1, i= j:degree+1, k = 1:q
        mass_line[j] += spline_val[i,k] * spline_val[i-j+1,k]*quad_w[k]
    end

	mass_line

end 

function ifft2d!( outval, solver, inval )

    nx, ny = solver.nx, solver.ny

    for i=1:nx
        for j=1:ny
            solver.array1d_y[j] = inval[i,j]
        end
        ifft!( solver.array1d_y )
        for j=1:ny
            inval[i,j] = solver.array1d_y[j]
        end
    end
    
    k=0

    for j=1:ny
        for i=1:nx
            solver.array1d_x(i) = inval(i,j)
        end
        ifft!(self.array1d_x)
        
        for i=1:nx
            k = k+1
            outval[k] = real( self.array1d_x[i] )
        end
    end

end
     

function fft2d( self, rho )

    k=0
    for j=1:ny
        for i=1:n_dofs(1)
            k = k+1
            solver.array1d_x[i] = rho[k]
        end
        fft!(solver.array1d_x)
        for i=1:n_dofs(1)
            solver.scratch[i,j] = solver.array1d_x[i]
        end
    end
    
       
    for i=1:nx
        for j=1:ny
            solver.array1d_y[j] = solver.scratch[i,j]
        end
        fft!(solver.array1d_y)
        for j=1:ny
            solver.scratch[i,j] = solver.array1d_y[j]
        end
    end

  
end 

# """
# Compute the L2 projection of a given function f on periodic splines of given degree
# """
# function l2projection_2d_fem(self, func, component, form)
# 
#     compute_fem_rhs(self, func, component, form, self.work )
# 
#     if form == 1 
#          solve(inverse_mass_1[component], self.work, coefs_dofs )
#     elseif form == 2
#          solve(inverse_mass_2[component], self.work, coefs_dofs )
#     end
# 
#     coefs_dofs
# 
# end 


function compute_e_from_rho_fft!( efield, solver,  rho )
  
  # Compute Fourier transform
  fft2d( solver, rho )
    
  # Apply inverse matrix of eigenvalues on mode
  for j=1:solver.nx, i=1:n_dofs[1]

      if ( i == 1 && j==1  )
         solver.scratch[i,j] = 0
      else
         eig_val = solver.eig_values_dtm1d_1[i] * solver.eig_values_mass_0_2[j] + solver.eig_values_mass_0_1[i] * solver.eig_values_dtm1d_2[j] 
         solver.scratch[i,j] = solver.scratch[i,j] / eig_val
      end
      
      solver.scratchx[i,j] = -solver.scratch[i,j]* solver.eig_values_d1[i]
      solver.scratchy[i,j] = -solver.scratch[i,j]* solver.eig_values_d2[j]
        
  end
  
    
  # Compute inverse Fourier transfrom
  

   ifft2d( solver, solver.scratchx, efield[1] )
   ifft2d( solver, solver.scratchy, efield[2] )
 
    
end

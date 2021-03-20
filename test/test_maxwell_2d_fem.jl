using Test

@testset "Maxwell 2D" begin

function evaluate_spline_2d( nx, ny, degs, dofs )

    deg1, deg2 = degs
    vals = zeros(nx * ny)
    a_in = zeros(ny) 
    a_out = zeros(ny)
    
    istart, iend = 1, nx
    for j=1:ny
       val = GEMPIC.eval_uniform_periodic_spline_curve(deg1, dofs[istart:iend])
       vals[istart:iend] .= val
       istart = iend+1
       iend = iend + nx
    end

    for i=1:nx
       for j=1:ny
          a_in[j] = vals[i+(j-1)*nx]
       end
       a_out .= GEMPIC.eval_uniform_periodic_spline_curve(deg2, a_in)
       for j=1:ny
          vals[i+(j-1)*nx] = a_out[j]
       end
    end

    @show size(vals)
    vals

end


  xmin, xmax = 0.0, 2π
  nx = 16
  ymin, ymax = 0.0, 2π
  ny = 32

  n1, n2 = nx, ny

  mesh = TwoDGrid( xmin, xmax, nx, ymin, ymax, ny)

  deg = 3
  delta_t = 0.01
  nsteps = 30
  
  maxwell = TwoDMaxwell(mesh, deg)

  x = zeros(nx+1,ny+1)
  y = zeros(nx+1,ny+1)
  efield = [ zeros(nx * ny) for _ in 1:3]
  bfield = deepcopy(efield)
  efield_ref = deepcopy(efield)
  bfield_ref = deepcopy(efield)

  x = LinRange(xmin, xmax, nx+1)[1:end-1] .* transpose(ones(ny))
  y = ones(nx) .* transpose(LinRange(xmin, xmax, ny+1)[1:end-1])

  w1 = sqrt(3)
  w2 = sqrt(3)
  
  nt = nx * ny
  rho = zeros( nt)
  rho_ref = zeros( nt) 
  time = 0.0

  b3(x, y) = - cos(x)*cos(y)*cos(sqrt(2)*time)
  e1(x, y) = cos(x)*sin(y)*sin(sqrt(2)*time)/sqrt(2)
  e2(x, y) = - sin(x)*cos(y)*sin(sqrt(2)*time)/sqrt(2)

  sin_k = (x, y) -> sin((x+y)-w1*time) 
  cos_k = (x, y) -> cos((x+y)-w1*time) 

  rho .= compute_rhs_from_function( maxwell, cos_k, 1, 0 )

  compute_e_from_rho!( efield, maxwell, rho )

  efield_val1 = evaluate_spline_2d( nx, ny, (deg-1,deg  ), efield[1])  
  efield_val2 = evaluate_spline_2d( nx, ny, (deg  ,deg-1), efield[2])
  efield_val3 = evaluate_spline_2d( nx, ny, (deg  ,deg  ), efield[3])
  
  efield_ref[1] .= vec(sin_k.(x, y) ./ 2)
  efield_ref[2] .= efield_ref[1]
  efield_ref[3] .= 0.0
 
  @test efield_val1 ≈ efield_ref[1] rtol=1e-4
  @test efield_val2 ≈ efield_ref[2] rtol=1e-4
  @test efield_val3 ≈ efield_ref[3]

  time = 0.0

  time = -0.5*delta_t

  bfield[1] .= l2projection( maxwell, e1, 1, 2)
  bfield[2] .= l2projection( maxwell, e2, 2, 2)
  bfield[3] .= l2projection( maxwell, b3, 3, 2)

  time = 0.0
  
  efield[1] .= l2projection( maxwell, e1, 1, 1)
  efield[2] .= l2projection( maxwell, e2, 2, 1)
  efield[3] .= l2projection( maxwell, b3, 3, 1)
  efield[3] .*= -1

  
  for istep = 1:nsteps
     compute_b_from_e!( bfield, maxwell, delta_t, efield )
     compute_e_from_b!( efield, maxwell, delta_t, bfield )
  end

  # Evaluate E and B at the grid points
  bfield_val1 = evaluate_spline_2d( nx, ny, (deg,deg-1), bfield[1])
  bfield_val2 = evaluate_spline_2d( nx, ny, (deg-1,deg), bfield[2])
  bfield_val3 = evaluate_spline_2d( nx, ny, (deg-1,deg-1), bfield[3])

  efield_val1 = evaluate_spline_2d( nx, ny, (deg-1,deg), efield[1])
  efield_val2 = evaluate_spline_2d( nx, ny, (deg,deg-1), efield[2])
  efield_val3 = evaluate_spline_2d( nx, ny, (deg,deg), efield[3])


  # Reference solutions
  time = (nsteps)*delta_t
  ind = 1
  for j = 1:n2
     for i = 1:n1
        efield_ref[1][ind] = e1(x[i,j], y[i,j])
        efield_ref[2][ind] = e2(x[i,j], y[i,j])
        efield_ref[3][ind] = -b3(x[i,j], y[i,j])
        ind = ind+1
     end
  end
  time = (nsteps-0.5)*delta_t
  ind = 1
  for j = 1:n2, i = 1:n1
      bfield_ref[1][ind] = e1(x[i,j], y[i,j])
      bfield_ref[2][ind] = e2(x[i,j], y[i,j])
      bfield_ref[3][ind] = b3(x[i,j], y[i,j])
      ind = ind+1
  end

#=

  @test efield_val1 ≈ efield_ref[1] rtol=1e-3
  @test efield_val2 ≈ efield_ref[2] rtol=1e-3
  @test efield_val3 ≈ efield_ref[3] rtol=1e-3

  @test bfield_val1 ≈ bfield_ref[1] rtol=1e-3
  @test bfield_val2 ≈ bfield_ref[2] rtol=1e-3
  @test bfield_val3 ≈ bfield_ref[3] rtol=1e-3


  l2projection( maxwell, cos_k, 1, 1, efield[1:nt] )
  
  error2 = maxwell%inner_product( efield(1:nt), efield(1:nt), 1, 1 ) - 2.0*pi^2
  println("Error in L2 norm squared: $error2")
  
  rho = compute_rhs_from_function( maxwell, sin_k, 1, 1 )

  ompute_e_from_j( maxwell, rho, 1, efield(1:nt) )
  evaluate_spline_2d( nc_eta, [deg-1,deg,deg], efield(1:nt), efield_val(1:nt))
  
  ind = 1
  for j = 1:nc_eta2, i = 1:nc_eta1)
      efield_ref[ind] = cos_k([x(i,j), y(i,j)]) - sin_k([x(i,j), y(i,j)])
      ind = ind+1
  end
  error5 = maximum(abs(efield_val(1:nt)-efield_ref(1:nt)))
  println("Error compute_e_from_j: $error")

  rho_ref = compute_rhs_from_function( maxwell, cos_k, 1, 0 )
  rho_ref = 2.0 .* rho_ref
  l2projection( maxwell, sin_k, 1, 1, efield(1:nt) )
  l2projection( maxwell, sin_k, 2, 1, efield(nt+1:nt*2) )
  l2projection( maxwell, sin_k, 3, 1, efield(nt*2+1:nt*3) )

  compute_rho_from_e!( rho, maxwell, efield )
  
  error6 =  maximum( abs.( rho .- rho_ref ) )
  println( " Error compute_rho_from_e: $error6 ")

=#

end 

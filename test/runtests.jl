using Test
using GEMPIC

"""
   test_maxwell_1d_fem()

   test 1D Maxwell spline finite element solver on a periodic grid

`L_x` domain dimensions and M is an integer.

```math
B_z(x,y,t) = \\cos(\\frac{2 M \\pi}{L_x} x)  \\cos(\\frac{2 M \\pi}{L_x} t)
```
```math
E_y(x,y,t) = \\sin(\\frac{2 M \\pi}{L_x} x)  \\sin(\\frac{2 M \\pi}{L_x} t)
```

"""
function test_maxwell_1d_fem()

  #sll_m_low_level_bsplines, only: sll_s_eval_uniform_periodic_spline_curve


  eta1_min = .0
  eta1_max = 2π
  nc_eta1  = 128

  Lx = eta1_max - eta1_min

  delta_eta1 = Lx / nc_eta1

  domain = [eta1_min, eta1_max]
  
  deg = 3

  maxwell_1d = Maxwell1DFEM(domain, nc_eta1, deg)

#=
  ! Allocate arrays
  SLL_CLEAR_ALLOCATE(ex(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ey(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(bz(1:nc_eta1),error)

  SLL_CLEAR_ALLOCATE(bz_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ex_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ey_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(rho(1:nc_eta1),error)

  SLL_CLEAR_ALLOCATE(sval(1:nc_eta1),error)


  ! Test Poisson
  !-------------
  ! Set exact solution
  do i = 1, nc_eta1
     xi = eta1_min + real(i-1,f64)*delta_eta1
     ex_exact(i) =   sin_k(xi)/(2.0_f64*mode*sll_p_pi/Lx)
  end do

  call maxwell_1d%compute_rhs_from_function( cos_k, deg, rho)

  call maxwell_1d%compute_e_from_rho( ex, rho ) 

  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_eval_uniform_periodic_spline_curve(deg-1, ex, sval)
  err_ex = maxval(sval-ex_exact)
  print*, 'error Poisson',  err_ex

  call sll_s_plot_two_fields_1d('ex',nc_eta1,sval,ex_exact,0,0.0_f64)
 
  ! Test Ampere
  !-------------
  ! Set time step
  dt = .5_f64 * delta_eta1
  ! Set exact solution
  do i = 1, nc_eta1
     xi = eta1_min + real(i-1,f64)*delta_eta1
     ex_exact(i) =   -cos_k(xi)*dt
  end do

  call maxwell_1d%compute_rhs_from_function(cos_k, deg-1, rho)
  ex = 0.0_f64
  call maxwell_1d%compute_E_from_j(dt*rho, 1, ex )

  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_eval_uniform_periodic_spline_curve(deg-1, ex, sval)
  err_ex2 = maxval(sval-ex_exact)
  print*, 'error Ampere',  err_ex2
  !call sll_plot_two_fields_1d('ex',nc_eta1,sval,ex_exact,0,0.0_f64)

  l2norm =  maxwell_1d%l2norm_squared(ex,deg-1)
  err_l2norm = l2norm - dt**2*sll_p_pi
  print*, 'error l2 norm', err_l2norm


  ! Test Maxwell on By and Ez 
  !--------------------------
  ! Set time stepping parameters
  time  = 0.0_f64
  dt    = .5_f64 * delta_eta1
  nstep = 10

  ! Compute initial fields 
  ex = 0.0_f64 ! 0-form -> splines of degree deg
  call maxwell_1d%L2projection(cos_k, deg-1, bz) ! 0-form -> splines of degree deg-1

  ! Time loop. Second order Strang splitting
  do istep = 1, nstep 
     call maxwell_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
     call maxwell_1d%compute_e_from_b(         dt, bz, ey)
     call maxwell_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
     
     time = time + dt

     do i = 1, nc_eta1
        xi = eta1_min + real(i-1,f64)*delta_eta1
        ey_exact(i) =   sin(mode*2*sll_p_pi*xi/Lx) * sin(mode*2*sll_p_pi*time/Lx)
        bz_exact(i) =   cos(mode*2*sll_p_pi*xi/Lx) * cos(mode*2*sll_p_pi*time/Lx)
     end do

     call sll_s_eval_uniform_periodic_spline_curve(deg, ey, sval)
     err_ey = maxval(sval-ey_exact)
     call sll_s_eval_uniform_periodic_spline_curve(deg-1, bz, sval)
     err_bz = maxval(sval-bz_exact)

     write(*,"(10x,' istep = ',I6)",advance="no") istep
     write(*,"(' time = ',g12.3,' sec')",advance="no") time
     write(*,"(' erreur L2 = ',2g15.5)") err_ey, err_bz

     call sll_s_plot_two_fields_1d('bz',nc_eta1,sval,bz_exact,istep,time)

  end do ! next time step

  tol = 1.0d-3
  if ((err_bz < tol) .and. (err_ey < tol) .and. (err_ex < tol) .and. (err_ex2 < tol) .and. (err_l2norm < tol)) then
     print*,'PASSED'
  endif

  call maxwell_1d%free()

  DEALLOCATE(ex)
  DEALLOCATE(ey)
  DEALLOCATE(bz)
  DEALLOCATE(bz_exact)
  DEALLOCATE(ex_exact)
  DEALLOCATE(ey_exact)
  DEALLOCATE(rho)

contains

  function cos_k(x)

    sll_real64             :: cos_k
    sll_real64, intent(in) :: x

    cos_k = cos(mode*2*sll_p_pi*x/Lx) 

  end function cos_k

  function sin_k(x)

    sll_real64             :: sin_k
    sll_real64, intent(in) :: x

    sin_k = sin(mode*2*sll_p_pi*x/Lx) 

  end function sin_k

=#

    true

end 

@test test_maxwell_1d_fem()

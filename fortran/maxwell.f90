module Maxwell

use zone
implicit none

integer, private :: i, j

real(kind=prec), private :: dex_dy, dey_dx
real(kind=prec), private :: dbz_dx, dbz_dy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( tm )

type( tm_mesh_fields ) :: tm

   !*** On utilise l'equation de Faraday sur un demi pas
   !*** de temps pour le calcul du champ magnetique  Bz 
   !*** a l'instant n puis n+1/2 apres deplacement des
   !*** particules

   do i=0,nx-1
   do j=0,ny-1
      dex_dy     = (tm%ex(i,j+1)-tm%ex(i,j)) / hy(j)
      dey_dx     = (tm%ey(i+1,j)-tm%ey(i,j)) / hx(i)
      tm%bz(i,j) = tm%bz(i,j) + 0.5 * dt * (dex_dy - dey_dx)
   end do
   end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere( tm )

type( tm_mesh_fields ) :: tm

   !*** Calcul du champ electrique E au temps n+1
   !*** sur les points internes du maillage
   !*** Ex aux points (i+1/2,j)
   !*** Ey aux points (i,j+1/2)

   do i=0,nx-1
   do j=1,ny-1
      dbz_dy = (tm%bz(i,j)-tm%bz(i,j-1)) / hhy(j)
      tm%ex(i,j) = tm%ex(i,j) + csq * dt * dbz_dy - dt * tm%jx(i,j)/e0
   end do
   end do

   do i=1,nx-1
   do j=0,ny-1
      dbz_dx = (tm%bz(i,j)-tm%bz(i-1,j)) / hhx(i)
      tm%ey(i,j) = tm%ey(i,j) - csq * dt * dbz_dx - dt * tm%jy(i,j)/e0
   end do
   end do

end subroutine ampere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine conditions_limites( tm, time )

type( tm_mesh_fields ) :: tm
real(kind=prec) :: a11,a12,a21,a22,b1,b2,dis
real(kind=prec) :: time, alpha, omega
integer :: mm=1 !parametre du cas entran

do i=0,nx-1
   dbz_dy = (tm%bz(i,0)-tm%bz(i,ny-1)) / hhy(0)
   tm%ex(i,0)  = tm%ex(i,0) + csq * dt * dbz_dy - dt * tm%jx(i,0)/e0
   tm%ex(i,ny) = tm%ex(i,0)
end do

do j=0,ny-1
   dbz_dx = (tm%bz(0,j)-tm%bz(nx-1,j)) / hhx(0)
   tm%ey(0,j)  = tm%ey(0,j) - csq * dt * dbz_dx - dt * tm%jy(0,j)/e0
   tm%ey(nx,j) = tm%ey(0,j)
end do

end subroutine conditions_limites

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine decalage( tm, tm1 )

type(tm_mesh_fields) :: tm, tm1

!*** Calcul des composantes des champs 
!*** sur les noeuds du maillage de rho
!*** par interpolation lineaire

do i=1,nx-1
   do j=1,ny-1
      tm1%ex(i,j) = ( hx(i)*tm%ex(i-1,j) + hx(i-1)*tm%ex(i,j) ) &
           & / (hx(i)+hx(i-1))
      tm1%ey(i,j) = ( hy(j)*tm%ey(i,j-1) + hy(j-1)*tm%ey(i,j) ) &
           & / (hy(j)+hy(j-1))
      tm1%bz(i,j) = ( ( hx(i)*tm%bz(i-1,j-1) + hx(i-1)*tm%bz(i,j-1) ) &
           & * hy(j) + ( hx(i)*tm%bz(i-1,j) + hx(i-1)*tm%bz(i,j) ) &
           & * hy(j-1) ) / ( (hx(i)+hx(i-1)) * (hy(j)+hy(j-1)) )
   end do
end do


do i=1,nx-1 !Sud et Nord
   tm1%ex(i,0) = ( hx(i)*tm%ex(i-1,0) + hx(i-1)*tm%ex(i,0) ) &
           & / (hx(i)+hx(i-1))
   tm1%ey(i,0) = ( hy(0)*tm%ey(i,ny-1) + hy(ny-1)*tm%ey(i,0) ) &
           & / (hy(0)+hy(ny-1))
   tm1%bz(i,0) = ( ( hx(i)*tm%bz(i-1,ny-1) + hx(i-1)*tm%bz(i,ny-1) ) &
           & * hy(0) + ( hx(i)*tm%bz(i-1,0) + hx(i-1)*tm%bz(i,0) ) &
           & * hy(ny-1) ) / ( (hx(i)+hx(i-1)) * (hy(0)+hy(ny-1)) )
   tm1%ex(i,ny) = tm1%ex(i,0) 
   tm1%ey(i,ny) = tm1%ey(i,0) 
   tm1%bz(i,ny) = tm1%bz(i,0) 
end do

do j=1,ny-1 !Ouest et Est
   tm1%ex(0,j) = ( hx(0)*tm%ex(nx-1,j) + hx(nx-1)*tm%ex(0,j) ) &
        & / (hx(0)+hx(nx-1))
   tm1%ey(0,j) = ( hy(j)*tm%ey(0,j-1) + hy(j-1)*tm%ey(0,j) ) &
        & / (hy(j)+hy(j-1))
   tm1%bz(0,j) = ( ( hx(0)*tm%bz(nx-1,j-1) + hx(nx-1)*tm%bz(0,j-1) ) &
        & * hy(j) + ( hx(0)*tm%bz(nx-1,j) + hx(nx-1)*tm%bz(0,j) ) &
        & * hy(j-1) ) / ( (hx(0)+hx(nx-1)) * (hy(j)+hy(j-1)) )
   tm1%ex(nx,j) = tm1%ex(0,j) 
   tm1%ey(nx,j) = tm1%ey(0,j) 
   tm1%bz(nx,j) = tm1%bz(0,j) 
end do

!Coins
tm1%ex(0,0) = ( hx(0)*tm%ex(nx-1,0) + hx(nx-1)*tm%ex(0,0) ) &
        & / (hx(0)+hx(nx-1))
tm1%ey(0,0) = ( hy(0)*tm%ey(0,ny-1) + hy(ny-1)*tm%ey(0,0) ) &
        & / (hy(0)+hy(ny-1))
tm1%bz(0,0) = ( ( hx(0)*tm%bz(nx-1,ny-1) + hx(nx-1)*tm%bz(0,ny-1) ) &
        & * hy(0) + ( hx(0)*tm%bz(nx-1,0) + hx(nx-1)*tm%bz(0,0) ) &
        & * hy(ny-1) ) / ( (hx(0)+hx(nx-1)) * (hy(0)+hy(ny-1)) )

tm1%ex(nx,0)  = tm1%ex(0,0) 
tm1%ex(nx,ny) = tm1%ex(0,0) 
tm1%ex(0,ny)  = tm1%ex(0,0) 

tm1%ey(nx,0)  = tm1%ey(0,0) 
tm1%ey(nx,ny) = tm1%ey(0,0) 
tm1%ey(0,ny)  = tm1%ey(0,0) 

tm1%bz(nx,0)  = tm1%bz(0,0) 
tm1%bz(nx,ny) = tm1%bz(0,0) 
tm1%bz(0,ny)  = tm1%bz(0,0) 

end subroutine decalage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Maxwell

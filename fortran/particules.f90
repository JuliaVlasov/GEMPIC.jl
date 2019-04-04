module particules

use zone
use quietstart

implicit none

integer, private :: ipart 
integer, private :: i, j

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dichotomie_x( a, b, R, eps ) 

! il faut D(a)<R<D(b), on cherche x tq R=D(x), resu dans a 

real(kind = prec) :: a, b, R, eps, x, D

do while ( .true. )
   x = (a+b)/2
   D = ( kx*x + alpha * sin(kx*x) ) / (2*pi)
   if ( D<R-eps ) then
      a = x
   else if ( D>R+eps ) then 
      b = x
   else
      a = x
      return
   endif
end do

end subroutine dichotomie_x


subroutine creapa( ele )

type (particle) :: ele
real(kind=prec) :: speed, theta, vth, n
real(kind=prec) :: a, b, eps, R
integer :: k

eps = 1.d-12

vth =  1.
nbpart = 100*nx*ny
n = 1.d0/nbpart

allocate(ele%pos(2,nbpart))
allocate(ele%cell(2,nbpart))
allocate(ele%vit(2,nbpart))
allocate(ele%epx(nbpart))
allocate(ele%epy(nbpart))
allocate(ele%bpz(nbpart))
allocate(ele%p(nbpart))

do k=0,nbpart-1

   speed = vth * sqrt(-2 * log( (k+0.5)*n ))

   theta = trinary_reversing( k ) * 2 * pi

   a = 0; b = dimx 
   R = bit_reversing( k )
   call dichotomie_x(a,b,R,eps) 
   ele%pos(1,k+1) = a
   ele%pos(2,k+1) = dimy * penta_reversing( k ) 

   i = 0
   do while (ele%pos(1,k+1) >= x(i)) 
      i=i+1
   enddo
   ele%cell(k+1,1) = i-1
   
   j = 0
   do while (ele%pos(2,k+1) >= y(j)) 
      j=j+1 
   enddo
   ele%cell(2,k+1) = j-1
   
   ele%vit(1,k+1) = speed * cos(theta)  !
   ele%vit(2,k+1) = speed * sin(theta)  !

   ele%p(k+1) = poids * n

enddo


end subroutine creapa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interpol_eb( tm1, ele )

type  (particle) :: ele
type(tm_mesh_fields) :: tm1
real(kind = prec) :: a1, a2, a3, a4
real(kind = prec) :: xp, yp, dum
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

do ipart=1,nbpart
   i = ele%cell(1,ipart)
   j = ele%cell(2,ipart)
   xp = ele%pos(1,ipart)
   yp = ele%pos(2,ipart)

   dum = 1./(hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum

   ele%epx(ipart) = a1 * tm1%ex(i,j) + a2 * tm1%ex(i+1,j) &
        & + a3 * tm1%ex(i+1,j+1) + a4 * tm1%ex(i,j+1) 
   ele%epy(ipart) = a1 * tm1%ey(i,j) + a2 * tm1%ey(i+1,j) &
        & + a3 * tm1%ey(i+1,j+1) + a4 * tm1%ey(i,j+1) 
   ele%bpz(ipart) =  a1 * tm1%bz(i,j) + a2 * tm1%bz(i+1,j) &
        & + a3 * tm1%bz(i+1,j+1) + a4 * tm1%bz(i,j+1) 
end do

end subroutine interpol_eb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_vitesse( ele )

type (particle) :: ele
real(kind=prec) :: dum, u2
real(kind=prec) :: tantheta, sintheta
real(kind=prec) :: gamma

do ipart = 1, nbpart

   !*** Changement de variable u = gamma*vit

   if( relativ ) then

      u2    =   ele%vit(1,ipart)*ele%vit(1,ipart)    &
              + ele%vit(2,ipart)*ele%vit(2,ipart)
      if ( u2 >= csq ) then 
         print*,'Erreur : u2 >= c2 dans le calcul de la vitesse'
         print*,'ipart = ',ipart,' vx = ',ele%vit(1,ipart),' vy = ',ele%vit(2,ipart)
         stop
      else
         gamma = 1./sqrt( 1. - u2/csq )
      endif

      ele%vit(1,ipart) = gamma*ele%vit(1,ipart)
      ele%vit(2,ipart) = gamma*ele%vit(2,ipart)

   else

      gamma=1.

   end if


   !*** Separation des effets electriques et magnetiques

   !*** On ajoute la moitie de l'effet champ electrique E

   dum = 0.5 * dt * q_sur_m
   ele%vit(1,ipart) = ele%vit(1,ipart) + dum*(ele%epx(ipart)+exext)
   ele%vit(2,ipart) = ele%vit(2,ipart) + dum*(ele%epy(ipart)+eyext)

   !*** Algorithme de Buneman pour les effets magnetiques
 
   tantheta = dum * (ele%bpz(ipart)+bzext) / gamma 
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   ele%vit(1,ipart) = ele%vit(1,ipart) + ele%vit(2,ipart)*tantheta
   ele%vit(2,ipart) = ele%vit(2,ipart) - ele%vit(1,ipart)*sintheta
   ele%vit(1,ipart) = ele%vit(1,ipart) + ele%vit(2,ipart)*tantheta

   !*** Autre moitie de l'effet du champ electrique E

   ele%vit(1,ipart) = ele%vit(1,ipart) + dum*(ele%epx(ipart)+exext)
   ele%vit(2,ipart) = ele%vit(2,ipart) + dum*(ele%epy(ipart)+eyext)

   !*** On repasse a la vitesse (changement de variable inverse)

   if( relativ ) then

      u2 =   ele%vit(1,ipart)*ele%vit(1,ipart)    &
           + ele%vit(2,ipart)*ele%vit(2,ipart)

      gamma = sqrt( 1. + u2/csq )
 
      ele%vit(1,ipart) = ele%vit(1,ipart) / gamma
      ele%vit(2,ipart) = ele%vit(2,ipart) / gamma

   end if

end do

end subroutine avancee_vitesse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_part( ele, coef )  !Avancee de coef * dt

type(particle) :: ele
real(kind=prec) :: coef

do ipart=1,nbpart     
   ele%pos(1,ipart) = ele%pos(1,ipart) + ele%vit(1,ipart)*dt*coef
   ele%pos(2,ipart) = ele%pos(2,ipart) + ele%vit(2,ipart)*dt*coef
enddo

!*** Mise a jour des "cells"

do ipart=1,nbpart
   i = 0
   do while (ele%pos(1,ipart) >= x(i) .and. ele%pos(1,ipart)<dimx) 
      i=i+1
   enddo
   if ( ele%pos(1,ipart) >= dimx ) i=nx+1
   ele%cell(1,ipart) = i-1
   j = 0
   do while (ele%pos(2,ipart) >= y(j) .and. ele%pos(2,ipart)<dimy) 
      j=j+1 
   enddo
   if ( ele%pos(2,ipart) >= dimy ) j=ny+1
   ele%cell(2,ipart) = j-1
end do

end subroutine avancee_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sortie_part( ele )

type(particle) :: ele

!*** Traitement de la sortie des particules

do ipart=1,nbpart
   if( ele%pos(1,ipart) >= dimx ) ele%pos(1,ipart) = ele%pos(1,ipart) - dimx
   if( ele%pos(2,ipart) >= dimy ) ele%pos(2,ipart) = ele%pos(2,ipart) - dimy
   if( ele%pos(1,ipart) < 0.d0 )  ele%pos(1,ipart) = ele%pos(1,ipart) + dimx
   if( ele%pos(2,ipart) < 0.d0 )  ele%pos(2,ipart) = ele%pos(2,ipart) + dimy
end do   

!*** Mise a jour des "cases"

do ipart=1,nbpart
   i = 0
   do while (ele%pos(1,ipart) >= x(i) .and. ele%pos(1,ipart)<=dimx) 
      i=i+1
   enddo
   ele%cell(1,ipart) = i-1
   j = 0
   do while (ele%pos(2,ipart) >= y(j) .and. ele%pos(2,ipart)<=dimy) 
      j=j+1 
   enddo
   ele%cell(2,ipart) = j-1
end do

end subroutine sortie_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcul_rho( ele, tm )

type(particle) :: ele
type(tm_mesh_fields) :: tm
real(kind=prec) :: a1, a2, a3, a4, dum, xp, yp
real(kind=prec) :: rho_total

tm%r0 = tm%r1   
tm%r1 = 0.d0    
                
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

do ipart=1,nbpart
   i = ele%cell(1,ipart)
   j = ele%cell(2,ipart)
   xp = ele%pos(1,ipart)
   yp = ele%pos(2,ipart)
   dum = ele%p(ipart) / (hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum
   tm%r1(i,j)     = tm%r1(i,j)     + a1/(hhx(i)*hhy(j)) !charge unite = 1
   tm%r1(i+1,j)   = tm%r1(i+1,j)   + a2/(hhx(i+1)*hhy(j)) 
   tm%r1(i+1,j+1) = tm%r1(i+1,j+1) + a3/(hhx(i+1)*hhy(j+1))
   tm%r1(i,j+1)   = tm%r1(i,j+1)   + a4/(hhx(i)*hhy(j+1))
end do

do i=0,nx
   tm%r1(i,0)  = tm%r1(i,0) + tm%r1(i,ny)
   tm%r1(i,ny) = tm%r1(i,0)
end do
do j=0,ny
   tm%r1(0,j)  = tm%r1(0,j) + tm%r1(nx,j)
   tm%r1(nx,j) = tm%r1(0,j)
end do

rho_total = 0.d0
do i=0,nx-1
   do j=0,ny-1
      rho_total = rho_total + tm%r1(i,j)*hhx(i)*hhy(j)
   enddo
enddo
print*,'rho total',rho_total 
! Neutralisation du milieu
tm%r1 = tm%r1 - rho_total/dimx/dimy

end subroutine calcul_rho

subroutine calcul_j_cic( ele, tm )

type(particle) :: ele
type(tm_mesh_fields) :: tm
real(kind=prec) :: a1, a2, a3, a4, dum, xp, yp
real(kind=prec), dimension(0:nx,0:ny) :: jx, jy

jx = 0.d0
jy = 0.d0

do ipart=1,nbpart
   i = ele%cell(1,ipart)
   j = ele%cell(2,ipart)
   xp = ele%pos(1,ipart)
   yp = ele%pos(2,ipart)
   dum = ele%p(ipart) / (hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum
   dum = ele%vit(1,ipart) / (hx(i)*hy(j)) !charge unite = 1
   jx(i,j)     = jx(i,j)     + a1*dum  
   jx(i+1,j)   = jx(i+1,j)   + a2*dum 
   jx(i+1,j+1) = jx(i+1,j+1) + a3*dum 
   jx(i,j+1)   = jx(i,j+1)   + a4*dum 
   dum = ele%vit(2,ipart) / (hx(i)*hy(j)) 
   jy(i,j)     = jy(i,j)     + a1*dum  
   jy(i+1,j)   = jy(i+1,j)   + a2*dum 
   jy(i+1,j+1) = jy(i+1,j+1) + a3*dum 
   jy(i,j+1)   = jy(i,j+1)   + a4*dum 
end do

do i=0,nx
   jx(i,0)  = jx(i,0) + jx(i,ny)
   jx(i,ny) = jx(i,0)
   jy(i,0)  = jy(i,0) + jy(i,ny)
   jy(i,ny) = jy(i,0)
end do
do j=0,ny
   jx(0,j)  = jx(0,j) + jx(nx,j)
   jx(nx,j) = jx(0,j)
   jy(0,j)  = jy(0,j) + jy(nx,j)
   jy(nx,j) = jy(0,j)
end do

do i=0,nx-1
do j=0,ny
   tm%jx(i,j) = 0.5 * (jx(i,j)+jx(i+1,j))
end do
end do

do i=0,nx
do j=0,ny-1
   tm%jy(i,j) = 0.5 * (jy(i,j)+jy(i,j+1))
end do
end do

end subroutine calcul_j_cic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module particules

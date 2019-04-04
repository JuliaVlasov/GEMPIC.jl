module villasenor

use zone

implicit none

integer, private :: i, j
integer, private :: ipart

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcul_j_villa( ele, tm )

type(particle) :: ele
type(tm_mesh_fields) :: tm
real(kind=prec) :: x1, y1, x2, y2
integer :: i1, j1, i2, j2

tm%jx = 0.; tm%jy = 0.


do ipart = 1, nbpart

   x1 = ele%pos(1,ipart) - dt * ele%vit(1,ipart)
   y1 = ele%pos(2,ipart) - dt * ele%vit(2,ipart)
   i = 0
   do while (x1 >= x(i) .and. x1<=dimx) 
      i=i+1
   end do
   i1 = i-1
   j = 0
   do while (y1 >= y(j) .and. y1<=dimy)
      j=j+1
   end do
   j1 = j-1
   
   if (i1>nx .or. j1>ny .or. i1<0 .or. j1<0) then
      print*,'erreur pour la cell de x(n)',i1,j1
      stop
   endif
   
   x2 = ele%pos(1,ipart)
   y2 = ele%pos(2,ipart)
   i2 = ele%cell(1,ipart)
   j2 = ele%cell(2,ipart)
   
   if ((i1.eq.i2).and.(j1.eq.j2)) then
     
      !--4 frontieres traversees
      
      call quatre_frontieres_period(i1,j1,x1,y1,x2,y2,tm,ele)
      
   else if ( (i2.eq.i1+1 .and. j2.eq.j1) &
        & .or. (i2.eq.i1-1 .and. j2.eq.j1) &
        & .or. (j2.eq.j1+1 .and. i2.eq.i1) &
        & .or. (j2.eq.j1-1 .and. i2.eq.i1) ) then
      
      !-- 7 frontieres_period traversees
      
      call sept_frontieres_period(i1,j1,i2,j2,x1,y1,x2,y2,tm,ele)
                
   else if ( (i2.eq.i1+1.and.j2.eq.j1+1) .or. &
        & (i2.eq.i1-1.and.j2.eq.j1-1) .or. &
        & (i2.eq.i1+1.and.j2.eq.j1-1) .or. &
        & (i2.eq.i1-1.and.j2.eq.j1+1) ) then   
      
      !-- 10 frontieres_period traversees
      
      call dix_frontieres_period(i1,j1,i2,j2,x1,y1,x2,y2,tm,ele)
                    
   else 
      stop 'calcul du courant non traite'
   endif
      
end do

!*** conditions aux limites periodiques pour Nord et Est sur J


do j = 0,ny-1
   tm%jy(nx,j) = tm%jy(0,j)
enddo

do i = 0,nx-1
   tm%jx(i,ny) = tm%jx(i,0)
enddo
   
end subroutine calcul_j_villa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine quatre_frontieres(i1,j1,xold,yold,xnew,ynew,tm,ele)

type(particle) :: ele
type(tm_mesh_fields) :: tm
integer :: i1, j1
real(kind = prec):: xold,yold,xnew,ynew
real(kind = prec) :: Jx1,Jx2,Jy1,Jy2,dum


if ( 0 <= i1 .and. i1 < nx .and. 0 <= j1 .and. j1 < ny ) then

   dum = ele%p(ipart)*(xnew-xold)/(hx(i1)*hy(j1)*dt)
   
   Jx1 = dum * (y(j1+1) - 0.5*(yold+ynew))   
   Jx2 = dum * (0.5*(yold+ynew) - y(j1))
   
   dum = ele%p(ipart)*(ynew-yold)/(hx(i1)*hy(j1)*dt)
   
   Jy1 = dum * (x(i1+1) - 0.5*(xold+xnew))
   Jy2 = dum * (0.5*(xold+xnew) - x(i1))
  
   
   tm%jx(i1,j1  ) = tm%jx(i1,j1  ) + Jx1 / hhy(j1)
   tm%jx(i1,j1+1) = tm%jx(i1,j1+1) + Jx2 / hhy(j1+1)


   tm%jy(i1  ,j1) = tm%jy(i1  ,j1) + Jy1 / hhx(i1)
   tm%jy(i1+1,j1) = tm%jy(i1+1,j1) + Jy2 / hhx(i1+1)

endif

end subroutine quatre_frontieres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sept_frontieres( i1, j1, i2, j2, x1, y1, x2, y2, tm, ele )

type(particle) :: ele
type(tm_mesh_fields) :: tm
integer :: i2, j2, i1, j1
real(kind = prec) :: x1, y1, x2, y2
real(kind = prec) :: xinter, yinter

if ((i2.eq.(i1+1)).and.(j2.eq.j1)) then

   xinter = x(i2)
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)
   call quatre_frontieres(i1,j1,x1,y1,xinter,yinter,tm,ele)
   call quatre_frontieres(i2,j2,xinter,yinter,x2,y2,tm,ele)

else if ((i2.eq.(i1-1)).and.(j2.eq.j1)) then

   xinter = x(i1)
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)
   call quatre_frontieres(i1,j1,x1,y1,xinter,yinter,tm,ele)
   call quatre_frontieres(i2,j2,xinter,yinter,x2,y2,tm,ele)

else if ((i2.eq.i1).and.(j2.eq.j1+1)) then

   yinter = y(j2)
   xinter = x1 + (x2-x1)/(y2-y1) * (yinter - y1)
   call quatre_frontieres(i1,j1,x1,y1,xinter,yinter,tm,ele)
   call quatre_frontieres(i2,j2,xinter,yinter,x2,y2,tm,ele)

else if ((i2.eq.i1).and.(j2.eq.j1-1)) then

   yinter = y(j1)
   xinter = x1 + (x2-x1)/(y2-y1) * (yinter - y1)
   call quatre_frontieres(i1,j1,x1,y1,xinter,yinter,tm,ele)
   call quatre_frontieres(i2,j2,xinter,yinter,x2,y2,tm,ele)

endif

end subroutine sept_frontieres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dix_frontieres( i1, j1, i2, j2, x1, y1, x2, y2, tm, ele ) 

type(particle) :: ele
type(tm_mesh_fields) :: tm
real(kind=prec) :: x2, y2, x1, y1, xinter, yinter
integer :: i2, j2, jinter
integer :: i1, j1

if (i2.eq.i1+1.and.j2.eq.j1+1) then
       
   xinter = x(i2)
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)
   j = 0
   do while ( yinter >= y(j) .and. yinter < dimy )
      j=j+1
   end do
   if (yinter >= dimy) j=ny+1
   jinter = j - 1
   
   if (jinter==j1) then             
      call quatre_frontieres(i1,j1,x1,y1,xinter,yinter,tm,ele)
      call sept_frontieres(i1+1,j1,i2,j2,xinter,yinter,x2,y2,tm,ele)
   else if (jinter==j1+1) then                 
      call sept_frontieres(i1,j1,i1,j1+1,x1,y1,xinter,yinter,tm,ele)
      call quatre_frontieres(i2,j2,xinter,yinter,x2,y2,tm,ele)
   endif

else if (i2.eq.i1-1.and.j2.eq.j1-1) then
   
   xinter = x(i1)   
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)       
   j = 0
   do while ( yinter >= y(j) .and. yinter < dimy )
      j=j+1
   end do 
   if (yinter >= dimy) j=ny+1
   jinter = j - 1

   if (jinter==j1) then             
      call quatre_frontieres(i1,j1,x1,y1,xinter,yinter,tm,ele)
      call sept_frontieres(i1-1,j1,i2,j2,xinter,yinter,x2,y2,tm,ele)
   else if (jinter==j1-1) then        
      call sept_frontieres(i1,j1,i1,j1-1,x1,y1,xinter,yinter,tm,ele)
      call quatre_frontieres(i2,j2,xinter,yinter,x2,y2,tm,ele)          
   endif
 
else if (i2.eq.i1+1.and.j2.eq.j1-1) then
   
   xinter = x(i2)   
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)            
   j = 0
   do while ( yinter > y(j) .and. yinter < dimy )
      j=j+1
   end do 
   if (yinter >= dimy) j=ny+1
   jinter = j - 1
   
   if (jinter==j1) then         
      call quatre_frontieres(i1,j1,x1,y1,xinter,yinter,tm,ele)
      call sept_frontieres(i1+1,j1,i2,j2,xinter,yinter,x2,y2,tm,ele)
   else if (jinter==j1-1) then                    
      call sept_frontieres(i1,j1,i1,j1-1,x1,y1,xinter,yinter,tm,ele)
      call quatre_frontieres(i2,j2,xinter,yinter,x2,y2,tm,ele)
   endif

else if (i2.eq.i1-1.and.j2.eq.j1+1) then
   
   xinter = x(i1)   
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)          
   j = 0
   do while ( yinter > y(j) .and. yinter < dimy )
      j=j+1
   end do
   if (yinter >= dimy) j=ny+1
   jinter = j - 1
   
   if (jinter==j1) then             
      call quatre_frontieres(i1,j1,x1,y1,xinter,yinter,tm,ele)
      call sept_frontieres(i1-1,j1,i2,j2,xinter,yinter,x2,y2,tm,ele)
   else if (jinter==j1+1) then            
      call sept_frontieres(i1,j1,i1,j1+1,x1,y1,xinter,yinter,tm,ele)
      call quatre_frontieres(i2,j2,xinter,yinter,x2,y2,tm,ele)           
   endif

endif

end subroutine dix_frontieres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      CL period                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine quatre_frontieres_period(i1,j1,xold,yold,xnew,ynew,tm,ele)

type(particle) :: ele
type(tm_mesh_fields) :: tm
integer :: i1,j1
real(kind = prec):: xold,yold,xnew,ynew
real(kind = prec) :: Jx1,Jx2,Jy1,Jy2,dum


dum = ele%p(ipart)*(xnew-xold)/(hx(modulo(i1,nx))*hy(modulo(j1,ny))*dt)
  
Jx1 = dum * (y(j1+1) - 0.5*(yold+ynew))   
Jx2 = dum * (0.5*(yold+ynew) - y(j1))
   
dum = ele%p(ipart)*(ynew-yold)/(hx(modulo(i1,nx))*hy(modulo(j1,ny))*dt)

Jy1 = dum * (x(i1+1) - 0.5*(xold+xnew))
Jy2 = dum * (0.5*(xold+xnew) - x(i1))


tm%jx(modulo(i1,nx),modulo(j1,ny)  ) = tm%jx(modulo(i1,nx),modulo(j1,ny)  ) &
     & + Jx1 / hhy(modulo(j1,ny))
tm%jx(modulo(i1,nx),modulo(j1+1,ny)) = tm%jx(modulo(i1,nx),modulo(j1+1,ny)) &
     & + Jx2 / hhy(modulo(j1+1,ny))


tm%jy(modulo(i1,nx)  ,modulo(j1,ny)) = tm%jy(modulo(i1,nx)  ,modulo(j1,ny)) &
     & + Jy1 / hhx(modulo(i1,nx))
tm%jy(modulo(i1+1,nx),modulo(j1,ny)) = tm%jy(modulo(i1+1,nx),modulo(j1,ny)) &
     & + Jy2 / hhx(modulo(i1+1,nx))

end subroutine quatre_frontieres_period

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sept_frontieres_period(i1, j1, i2, j2, x1, y1, x2, y2, tm, ele)

type(particle) :: ele
type(tm_mesh_fields) :: tm
integer :: i1, j1, i2, j2
real(kind = prec) :: x1, y1, x2, y2
real(kind = prec) :: xinter,yinter

if ((i2.eq.(i1+1)).and.(j2.eq.j1)) then

   xinter = x(i2)
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)
   call quatre_frontieres_period(i1,j1,x1,y1,xinter,yinter,tm,ele)
   call quatre_frontieres_period(i2,j2,xinter,yinter,x2,y2,tm,ele)

else if ((i2.eq.(i1-1)).and.(j2.eq.j1)) then

   xinter = x(i1)
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)
   call quatre_frontieres_period(i1,j1,x1,y1,xinter,yinter,tm,ele)
   call quatre_frontieres_period(i2,j2,xinter,yinter,x2,y2,tm,ele)

else if ((i2.eq.i1).and.(j2.eq.j1+1)) then

   yinter = y(j2)
   xinter = x1 + (x2-x1)/(y2-y1) * (yinter - y1)  
   call quatre_frontieres_period(i1,j1,x1,y1,xinter,yinter,tm,ele)
   call quatre_frontieres_period(i2,j2,xinter,yinter,x2,y2,tm,ele)

else if ((i2.eq.i1).and.(j2.eq.j1-1)) then

   yinter = y(j1)
   xinter = x1 + (x2-x1)/(y2-y1) * (yinter - y1)
   call quatre_frontieres_period(i1,j1,x1,y1,xinter,yinter,tm,ele)
   call quatre_frontieres_period(i2,j2,xinter,yinter,x2,y2,tm,ele)

endif

end subroutine sept_frontieres_period

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dix_frontieres_period( i1, j1, i2, j2, x1, y1, x2, y2, tm, ele )

type(particle) :: ele
type(tm_mesh_fields) :: tm
real(kind=prec) :: x2, y2, x1, y1, xinter, yinter
integer :: i2, j2, jinter
integer :: i1, j1

if (i2.eq.i1+1.and.j2.eq.j1+1) then
       
   xinter = x(i2)
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)
   j = 0
   do while ( yinter >= y(j) .and. yinter < dimy )
      j=j+1
   end do
   if (yinter >= dimy) j=ny+1
   jinter = j - 1
   
   if (jinter==j1) then             
      call quatre_frontieres_period(i1,j1,x1,y1,xinter,yinter,tm,ele)
      call sept_frontieres_period(i1+1,j1,i2,j2,xinter,yinter,x2,y2,tm,ele)
   else if (jinter==j1+1) then                 
      call sept_frontieres_period(i1,j1,i1,j1+1,x1,y1,xinter,yinter,tm,ele)
      call quatre_frontieres_period(i2,j2,xinter,yinter,x2,y2,tm,ele)
   endif

else if (i2.eq.i1-1.and.j2.eq.j1-1) then
   
   xinter = x(i1)   
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)       
   j = 0
   do while ( yinter >= y(j) .and. yinter < dimy )
      j=j+1
   end do 
   if (yinter >= dimy) j=ny+1
   jinter = j - 1

   if (jinter==j1) then             
      call quatre_frontieres_period(i1,j1,x1,y1,xinter,yinter,tm,ele)
      call sept_frontieres_period(i1-1,j1,i2,j2,xinter,yinter,x2,y2,tm,ele)
   else if (jinter==j1-1) then        
      call sept_frontieres_period(i1,j1,i1,j1-1,x1,y1,xinter,yinter,tm,ele)
      call quatre_frontieres_period(i2,j2,xinter,yinter,x2,y2,tm,ele)          
   endif
 
else if (i2.eq.i1+1.and.j2.eq.j1-1) then
   
   xinter = x(i2)   
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)            
   j = 0
   do while ( yinter > y(j) .and. yinter < dimy )
      j=j+1
   end do 
   if (yinter >= dimy) j=ny+1
   jinter = j - 1
   
   if (jinter==j1) then         
      call quatre_frontieres_period(i1,j1,x1,y1,xinter,yinter,tm,ele)
      call sept_frontieres_period(i1+1,j1,i2,j2,xinter,yinter,x2,y2,tm,ele)
   else if (jinter==j1-1) then                    
      call sept_frontieres_period(i1,j1,i1,j1-1,x1,y1,xinter,yinter,tm,ele)
      call quatre_frontieres_period(i2,j2,xinter,yinter,x2,y2,tm,ele)
   endif

else if (i2.eq.i1-1.and.j2.eq.j1+1) then
   
   xinter = x(i1)   
   yinter = y1 + (y2-y1)/(x2-x1) * (xinter - x1)          
   j = 0
   do while ( yinter > y(j) .and. yinter < dimy )
      j=j+1
   end do
   if (yinter >= dimy) j=ny+1
   jinter = j - 1
   
   if (jinter==j1) then             
      call quatre_frontieres_period(i1,j1,x1,y1,xinter,yinter,tm,ele)
      call sept_frontieres_period(i1-1,j1,i2,j2,xinter,yinter,x2,y2,tm,ele)
   else if (jinter==j1+1) then            
      call sept_frontieres_period(i1,j1,i1,j1+1,x1,y1,xinter,yinter,tm,ele)
      call quatre_frontieres_period(i2,j2,xinter,yinter,x2,y2,tm,ele)           
   endif

endif

end subroutine dix_frontieres_period

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module villasenor

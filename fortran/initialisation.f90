module initialisation

use zone

implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initialisation des champs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init( tm )

type (tm_mesh_fields) :: tm

real(kind=prec) :: aux1, aux2
integer :: i, j

tm%ex = 0.d0; tm%ey = 0.d0; tm%bz = 0.d0
tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0
do i=0,nx-1
   aux1 = alpha/kx * sin(kx*x(i))
   aux2 = alpha * cos(kx*x(i))
   do j=0,ny
      tm%ex(i,j) = aux1
      tm%r1(i,j) = aux2
   enddo
enddo
      
end subroutine init

end module initialisation

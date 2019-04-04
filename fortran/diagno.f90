module diagno

use zone

implicit none


contains

subroutine modeE( tm, iplot, time )

type(tm_mesh_fields), intent(in) :: tm
integer             , intent(in) :: iplot
real(kind=prec)     , intent(in) :: time

real(kind=prec) :: aux
integer         :: i, j

aux =0.d0
do i=0,nx-1
   do j=0,ny-1
      aux = aux + tm%ex(i,j)*tm%ex(i,j)*hx(i)*hhy(j)
   end do
end do
aux = 0.5*log(aux)

open(34,file='modeE.dat',position="append")
if (iplot==1) rewind(34)
write(34,*) time, aux
close(34)

end subroutine modeE

end module diagno

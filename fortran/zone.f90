module zone

integer, parameter :: prec=8

type tm_mesh_fields
   real(kind=prec), dimension(:,:), pointer :: ex, ey
   real(kind=prec), dimension(:,:), pointer :: bz
   real(kind=prec), dimension(:,:), pointer :: r0, r1
   real(kind=prec), dimension(:,:), pointer :: jx, jy
end type tm_mesh_fields

type particle
   real(kind=prec)   , pointer :: pos(:,:)
   integer           , pointer :: cell(:,:)
   real(kind=prec)   , pointer :: vit(:,:)
   real(kind=prec)   , pointer :: epx(:)
   real(kind=prec)   , pointer :: epy(:)
   real(kind=prec)   , pointer :: bpz(:)
   real(kind=prec)   , pointer :: p(:)
end type particle

logical :: relativ 

real(kind=prec) :: pi 

character(len=6) :: jname


integer :: nx, ny
integer :: nstep, nstepmax
integer :: icrea, idiag
integer :: nbpart

integer, private :: i, j

real(kind=prec) :: dt, alpha, kx, ky, c, csq, e0
real(kind=prec), private :: dx, dy
real(kind=prec), dimension(:), allocatable :: x, y
real(kind=prec), dimension(:), allocatable :: hx, hy    ! les h_i+1/2
real(kind=prec), dimension(:), allocatable :: hhx, hhy  ! les h_i
real(kind=prec) :: dimx, dimy, dimx1, dimy1, dimx2, dimy2
real(kind=prec) :: cfl
real(kind=prec) :: tfinal

real(kind=prec) :: exext, eyext, bzext

real(kind=prec) :: charge, masse, poids
real(kind=prec) :: q_sur_m 


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readin( )

implicit none

pi = 4. * atan(1.)

alpha = 0.1
kx = 0.5
ky = 0.

dimx      = 2*pi/kx
dimy      = 1.             !dimensions du domaine 
nx        = 120            !nombre de pts suivant x
ny        = 10             !nombre de pts suivant y
cfl       = 0.9            !nombre de Courant-Friedrich-Levy
tfinal    = 20.            !temps final
nstepmax  = 20000000    !nbre d'iterations maxi
jname     = 'jcico1'       !methode de calcul des courants
icrea     = 1              !frequence d'emission des particules
idiag     = 10             !frequence des sorties graphiques
exext     = 0.             !champ electrique exterieur suivant x
eyext     = 0.             !champ electriaue exterieur suivant y
bzext     = 0.             !champ magnetique exterieur
charge    = 1.d0           !charge d'une macro particule
masse     = 1.d0           !masse d'une macro particule
c         = 8.             !vitesse de la lumiere
e0        = 1.             !permittivite du vide
relativ   = .false.

csq = c*c
q_sur_m = charge / masse
poids = charge

poids = dimx * dimy ! car int(f0) = dimx*dimy

!Creation du maillage

allocate(x(-1:nx+1))  !0:nx))
allocate(y(-1:ny+1))  !0:ny))
allocate(hx(-1:nx))
allocate(hy(-1:ny))
allocate(hhx(0:nx))
allocate(hhy(0:ny))

dx = dimx / nx
dy = dimy / ny

x(0) = 0.
y(0) = 0.

do i=1,nx
   x(i) = i*dx !(i*dx) *(i*dx+1)/(1+dimx)
enddo
do j=1,ny
   y(j) = j*dy ! (j*dy) *(j*dy+1)/(1+dimy)
enddo

do i=0,nx-1
   hx(i) = x(i+1)-x(i)
end do
do j=0,ny-1
   hy(j) = y(j+1)-y(j)
end do
hx(nx) = hx(0)  ! CL periodiques
hx(-1) = hx(nx-1)
hy(ny) = hy(0)
hy(-1) = hy(ny-1)

x(-1)   = x(0) - hx(nx-1)  !points utiles pour le cas period
x(nx+1) = x(nx) + hx(0)
y(-1)   = y(0) - hy(ny-1)
y(ny+1) = y(ny) + hy(0)

hhx(0) =  0.5 * ( hx(0) + hx(nx-1) )  !0.5 * hx(0)
hhx(nx) =  0.5 * ( hx(0) + hx(nx-1) )   !0.5 * hx(nx-1)
do i=1,nx-1
   hhx(i) = 0.5 * ( hx(i) + hx(i-1) )
enddo
hhy(0) = 0.5 * ( hy(0) + hy(ny-1) )   !0.5 * hy(0)
hhy(ny) = 0.5 * ( hy(0) + hy(ny-1) )   !0.5 * hy(ny-1)
do j=1,ny-1
   hhy(j) = 0.5 * ( hy(j) + hy(j-1) )
enddo

dx = hx(0)  !on calcule le plus petit pas
do i=1,nx-1
   if (hx(i)<dx)  dx = hx(i)
end do

dy = hy(0)
do j=1,ny-1
   if (hy(j)<dy)  dy = hy(j)
end do

dt    = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c

nstep = floor(tfinal/dt)

write(*,*) " cfl = ", cfl
write(*,*) " dx = ", dx, " dy = ", dy, " dt = ", dt
if( nstep > nstepmax ) nstep = nstepmax
write(*,*) " Nombre d'iteration nstep = ", nstep
write(*,*)

end subroutine readin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zone

!Ce code Fortran modelise le mouvement d'une planete autour du soleil. GM=1, M>>m

!Execution:
!gfortran equadiff.f90 -o equadiff.x

!Gnuplot:
!plot 'rk4.out' using 2:3
!set xlabel 'x'
!set ylabel 'y'
!replot

!_________________________

!X(1)=x(t)
!X(2)=y(t)
!X(3)=vx(t)
!X(4)=vy(t)

!Xderiv(1)=vx(t)
!Xderiv(2)=vy(t)
!Xderiv(3)=-x/(r^3)
!Xderiv(4)=-y/(r^3)

!conditions initiales annoncees dans l'enonce // 
!x0=0 // x0=1
!y0=1 // y0=0
!vx0=-0.5 // vx0=0
!vy0=0 // vy0=-0.5

program  meca_planet
implicit none
real (8) :: dt ,t ,ener
integer  :: i
real(8),  dimension (4) :: X, Xderiv
external  :: euler, deriv_planet, rk4

t=0.
dt=0.01

!Initialization
X=(/1.,0.,0.,-0.5/)

open(11,file='rk4.out')

do i=1,100000!tmax=0.01*100000=1000
	t=t+dt
	call  rk4(t,X,dt,4,deriv_planet)
	ener =0.5*(X(3)**2+X(4)**2)-1/sqrt(X(1)**2+X(2)**2)
	if (mod(nint(t/dt),10).eq.0) then
		write(11,*) t, X!, ener
	endif
enddo
close (11)
end program meca_planet

subroutine deriv_planet(t,X,Xderiv,n)
implicit none
integer, intent(in) :: n
real(8), intent (in) :: t
real(8) :: xnew, ynew, radius
real(8), dimension(n), intent(in) :: X
real(8), dimension(n), intent(out) :: Xderiv
if (n.ne.4) write (*,*) 'WARNING: dimension de n incorrecte, devrait etre 4'
radius=sqrt(X(1)**2+X(2)**2)
Xderiv(1)=X(3)
Xderiv(2)=X(4)
Xderiv(3)=-X(1)/radius**3
Xderiv(4)=-X(2)/radius**3
end subroutine deriv_planet

subroutine rk4(t,X,dt,n,deriv)
!Runge-Kutta du 4eme ordre
implicit none
integer, intent(in) :: n
real(8), intent(in) :: t, dt
real(8), dimension(n), intent (inout) :: X
real(8) :: ddt
real(8), dimension(n) :: Xp, k1, k2, k3, k4
ddt = 0.5*dt
call deriv(t,X,k1,n); Xp=X+ddt*k1
call deriv(t+ddt,Xp,k2,n); Xp=X+ddt*k2
call deriv(t+ddt,Xp,k3,n); Xp=X+dt*k3
call deriv(t+dt,Xp,k4,n); X = X + dt*(k1+2.0*k2 + 2.0*k3+k4)/6.0
end subroutine rk4

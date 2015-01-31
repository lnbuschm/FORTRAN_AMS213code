!  ********** HOMEWORK #5 --------- AMS 213 *************** !
!  Luke Buschmann
program hw5
use LinAl
use DiffEq
implicit none
real :: deriv2
character (len=100) :: filename
integer :: i,j,k
real :: x=1.,y,guess=1.,dt,t0,tf
!external deriv,deriv2
	!integer, parameter :: order = 3

! Study the stability of Euler's modified algorithm. Is it always stable, conditionally stable, or always unsable? 
!   If it is conditionally stable, what is the condition on the step for stability?
! Ans:  Euler's modified algorithm can be numerically unstable (it is conditionally stable)
! If the Euler method is applied to the linear equation y' = k y, then the numerical solution is unstable 
!   if the product hk is outside the region  \{ z \in \mathbf{C} \mid |z+1| \le 1 \},
!  z+1 <= 1 
! Conditionally stable for step sizes 0 < dt < 2 . ... 

! Study the stability of the Adams-Bashforth 2nd order scheme. Is it always stable, conditionaly stable, or always unsable?
!   If it is conditionally stable, what is the condition on the step for stability?
!  Conditionally stable for 0 < dt < 1

! Implement both schemes numerically. This time, the scheme must be able to deal with any number of coupled ODEs. 
!   As in the previous homework: create a timestepper routine, and a driver routine. Make sure the timestepper routine
!   takes in as argument the RHS function/routine that returns the derivatives. Bonus point for elegant use of derived data types.

! Compare the solutions obtained with the same timestep with both schemes, with the exact solution, for the problem
!    f'' = -f' -f, with initial conditions f(0) = 2, f'(0) = 1. Test different timesteps.
!  solve ordinary differential equation y''(t)=-y'(t)-y(t), y(0)=2, y'(0)=1
!     y1 = f,  y2 = f',          y1' = y2   ,   y2' = -y2 - y1
! Algebraic Solution:  y(t) = 2.0/3.0 * exp(-t/2.) * (2.0 * sqrt(3.) * sin(sqrt(3.) * t / 2.) + 3. * cos( sqrt(3.) * t / 2.))
!  y(t) = 0.666667 * 2.71828^(-0.5 * t) * (3. * cos(0.866025 * t) + 3.4641 * sin(0.866025 * t))
!         0.666667*  2.71828**(-0.5* x)* (3.0* cos(0.866025*x)+3.4641* sin(0.866025* x)) 
! Gnuplot:  
!   plot [0:10] 2.0/3.0 * exp(-x/2.) * (2.0 * sqrt(3.) * sin(sqrt(3.) * x / 2.) + 3. * cos( sqrt(3.) * x / 2.))
!y(t) = 
!  plot 'C:\cygwin64\home\Luke\eulermodified.dat' 
! plot [0:4] 'C:\cygwin64\home\Luke\eulermodified.dat',  2.0/3.0 * exp(-x/2.) * (2.0 * sqrt(3.) * sin(sqrt(3.) * x / 2.) + 3. * cos( sqrt(3.) * x / 2.))
! plot [0:100] 'C:\cygwin64\home\Luke\eulermodified.dat',  2.0/3.0 * exp(-x/2.) * (2.0 * sqrt(3.) * sin(sqrt(3.) * x / 2.) + 3. * cos( sqrt(3.) * x / 2.))
!  plot [0:100] 'C:\cygwin64\home\Luke\adamsbashforth2.dat',  2.0/3.0 * exp(-x/2.) * (2.0 * sqrt(3.) * sin(sqrt(3.) * x / 2.) + 3. * cos( sqrt(3.) * x / 2.))

write(*,*) '******** HOMEWORK #5  --  AMS213 ******** '

dt = 0.5
!  Homework Note:
!   With dt > 1, the adams bashforth 2nd order algorithm is unstable. The points diverge from the ODE.
!   With dt > 2, both the euler modified and adams bashforth 2nd order become unstable. 
t0 = 0.0
tf = 100.0

!write(*,*)  'Enter Timestep: '
!read(*,*) dt
!write(*,*)  'Enter Initial Time: '
!read(*,*) t0
!write(*,*)  'Enter Final Time: '
!read(*,*) tf
	
!call ODEDRIVER(deriv,t0,tf,dt,'eulerexplicit.dat','EulerExplicit')

call ODEDRIVER(deriv,t0,tf,dt,'eulermodified.dat','EulerModified')

call ODEDRIVER(deriv,t0,tf,dt,'adamsbashforth2.dat','AdamsBashforth2')

end program hw5
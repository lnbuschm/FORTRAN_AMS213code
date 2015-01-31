!  ********** HOMEWORK #4 --------- AMS 213 *************** !
!  Luke Buschmann
program hw4
use LinAl
use DiffEq
implicit none
integer, parameter :: order = 1
real :: deriv1,deriv2,dt,t0,tf
REAL, POINTER :: deriv(:, :)
!real, dimension(:), allocatable :: deriv
character (len=100) :: filename
integer :: i,j,k
external deriv1,deriv2

!    Create a routine that advances the solution to an ODE of the kind df/dt = G(f,t) for one timestep using Euler's method. 
!     The input and outputs to the routine must be a derived data type that contains (1) the time and (2) the value of the
!     function f at that time. The name of the function that contains G(f,t) must be passed as an argument to the routine.

!    Create a routine that advances the solution to an ODE of the kind df/dt = G(f,t) for one timestep using Euler's modified 
!     midpoint method. The input and outputs to the routine must be a derived data type that contains (1) the time and (2)
!     the value of the function f at that time. The name of the function that contains G(f,t) must be passed as an argument
!     to the routine.

!    Create a driver program that integrates the ODE above using one method, and then the other. The program must take in as
!     arguments the initial time, the initial value of f, the final time, and the timestep. The user must know how to modify 
!     the routine that calculates G(f,t).

!    In Section, you will be given an ODE, an initial time, an initial f, and a final time, and you will have to produce 2 files
!  that contain the solutions obtained in both cases. You will be free to adjust the timestep as needed for a decent solution,
!  but you will have to use the same timestep in both cases. If you wish, you may put the routines in a module.

!allocate(deriv(order,order))
!deriv(1)=>deriv1

write(*,*) '******** HOMEWORK #4  --  AMS213 ******** '

!write(*,*)  'Enter Timestep: '
!read(*,*) dt
!write(*,*)  'Enter Initial Time: '
!read(*,*) t0
!write(*,*)  'Enter Final Time: '
!read(*,*) tf	
	
!call ODEDRIVER(deriv1,t0,tf,dt,'eulerexplicitfit.dat','EulerExplicit')
!call ODEDRIVER(deriv1,t0,tf,dt,'eulermodifiedfit.dat','EulerModified')


! In-class quiz
! Using y' = y - y^3,  f(0)=.1,  t0=0, tf=10, dt=0.1
call ODEDRIVER(deriv1,0.0,10.0,0.1,'eulerexplicitfit.dat','EulerExplicit')
call ODEDRIVER(deriv1,0.0,10.0,0.1,'eulermodifiedfit.dat','EulerModified')

end program hw4

! y' = y - y^3
function deriv1(x,y)
implicit none
real :: deriv1,x,y
deriv1 = y - y**3
end function deriv1

! y' = y + yx
function deriv2(x,y)
implicit none
real :: deriv2,x,y
deriv2 = y + x * y
end function deriv2


module DiffEq
implicit none
integer :: nsize,msize
integer,private :: i,j,k
type state ! state of ode
	real :: t
	real, dimension(:), allocatable :: x
	real, dimension(:,:), allocatable :: y
	integer :: step
endtype state
integer, parameter :: order = 2   ! Order of ODE in deriv() - input parameter
real, dimension(:), allocatable :: y0
REAL, PARAMETER :: Pi = 3.1415927
contains

! Subroutine CRANK_NICHOLSON_2D_INIT is the initialization routine for the 
!   Crank-Nicholson solver for the 2D diffusion equation. It creates the A matrix
!   and performs a cholesky decomposition for use in the solve routine.
! 
! ** Input parameters **
! x : input x vector (evenly spaced)
! y : input y vector (evenly spaced)
! bigI : input mesh I value
! bigJ : input mesh J value
! dt : timestep
! D : diffusion coefficient
! fg1,fg2,fg3,fg4 : dirichlet boundary conditions
! f3a,f3b : advection boundary functions
! A : out matrix 'U' cholesky decomposed matrix of LHS (n+1 timestep)
SUBROUTINE CRANK_NICHOLSON_2D_INIT(x,y,bigI,bigJ,dt,D,fg1,fg2,fg3,fg4,A)
	use LinAl
	IMPLICIT NONE
	integer :: i,j,k,m,n,bigI,bigJ,INFO
	real :: dy,dx,dxf,dyf,dt,D
	real, dimension(0:bigI+1) :: x
	real, dimension(0:bigJ+1) :: y
	real, dimension(1:bigI*bigJ, 1:bigI*bigJ) :: A
	real, dimension(bigI*bigJ,1) :: B
	real, dimension(0:bigI+1,0:bigJ+1,0:2) :: u
	real, dimension(bigI*bigJ) :: w
	character (len=100) :: filename = 'hw6debugA.dat'   ! Output data filename
	character (len=100) :: filename2 = 'hw6debugB.dat'   ! Output data filename
	real :: fg1,fg2,fg3,fg4
	external fg1,fg2,fg3,fg4
	
	A = 0.0
	B = 0.0
	
	! Calculate dy, dx (assume x and y vectors are evenly spaced)
	dy = y(1) - y(0)
	dx = x(1) - x(0)
	dyf = 1. / (dy ** 2)
	dxf = 1. / (dx ** 2)
	
	! Create A and RHS vector B
	do j=1,bigJ 
		do i=1,bigI  
			! calculate n index
			k = i + (j-1) * bigI
			
			A(k,k) = 1.0 + dt * D * (dxf + dyf) 
			
			! Check for indexing to boundary conditions
			if (i .lt. bigI) then
				A(k,k+1) = -dxf * dt * D / 2.0
			endif
				
			if (i .gt. 1) then
				A(k,k-1) =  -dxf * dt * D / 2.0
			endif
			
			if (j .lt. bigJ) then
				A(k,k+bigI) = -dyf * dt * D / 2.0
			endif
			
			if (j .gt. 1) then
				A(k,k-bigI) = -dyf * dt * D / 2.0
			endif

		enddo
	enddo

	!  For Debug: Output matrix A and B
	!call writematrix(A,filename)
	!call writematrix(B,filename2)
	
!	Call LAPACK routine to solve cholesky decomposition
	call SPOTRF( 'U', size(a,1), A, size(a,1), INFO )
	if (INFO .ne. 0) write(*,*) 'ERROR CALCULATING CHOLESKY DECOMPOSITION.... FAILED',INFO
 
END SUBROUTINE CRANK_NICHOLSON_2D_INIT


! Subroutine CRANK_NICHOLSON_2D_SOLVER is a Crank-Nicholson solver for the 2D diffusion equation
! 
! ** Input parameters **
! x : input x vector (evenly spaced)
! y : input y vector (evenly spaced)
! bigI : input mesh I value
! bigJ : input mesh J value
! dt : timestep
! D : diffusion coefficient
! fg1,fg2,fg3,fg4 : dirichlet boundary conditions
! f3a,f3b : advection boundary functions
! A : 'U' cholesky decomposed matrix of LHS (n+1 timestep)
! u : output solution mesh matrix
SUBROUTINE CRANK_NICHOLSON_2D_SOLVER(x,y,bigI,bigJ,dt,D,fg1,fg2,fg3,fg4,A,u)
	use LinAl
	IMPLICIT NONE
	integer :: i,j,k,m,n,bigI,bigJ,INFO
	real :: dy,dx,dxf,dyf,dt,D
	real, dimension(0:bigI+1) :: x
	real, dimension(0:bigJ+1) :: y
	real, dimension(1:bigI*bigJ, 1:bigI*bigJ) :: A
	real, dimension(bigI*bigJ,1) :: R ! B
	real, dimension(0:bigI+1,0:bigJ+1,0:2) :: u
	real, dimension(bigI*bigJ) :: w
	character (len=100) :: filename = 'hw6debugA.dat'   ! Output data filename
	character (len=100) :: filename2 = 'hw6debugB.dat'   ! Output data filename
	real :: fg1,fg2,fg3,fg4
	external fg1,fg2,fg3,fg4
	
	R = 0.0
	
	! Calculate dy, dx (assume x and y vectors are evenly spaced)
	dy = y(1) - y(0)
	dx = x(1) - x(0)
	dyf = 1. / (dy ** 2)
	dxf = 1. / (dx ** 2)
	
	do i=1,bigI
		do j=1,bigJ
			k = i + (j-1) * bigI
			
			R(k,1) = u(i,j,0) * (1 - (dt * D * (dxf + dyf))) !  u(i,j,0) - u(i,j,0) * dt * D * (dxf + dyf)

			! Check for indexing to boundary conditions
			if (i .eq. bigI) then ! i=bigI
				R(k,1) = R(k,1) + dt * D * dxf * fg4(y(j))
			else  ! i < bigI
				R(k,1) = R(k,1) + dt * D * dxf / 2.0 * u(i+1,j,0)
			endif
			
			if (i .eq. 1) then ! i=1
				R(k,1) = R(k,1) + dt * D * dxf  * fg3(y(j))
			else ! i > 1
				R(k,1) = R(k,1) + dt * D * dxf / 2.0 * u(i-1,j,0)
			endif
			
			
			if (j .eq. bigJ) then  ! j=bigJ
				R(k,1) = R(k,1) + dt * D * dyf  * fg2(x(i))
			else  ! j < bigJ
				R(k,1) = R(k,1) + dt * D * dyf / 2.0 * u(i,j+1,0)
			endif
			
			if (j .eq. 1) then ! j=1
				R(k,1) = R(k,1) + dt * D * dyf  * fg1(x(i))
			else ! j > 1
				R(k,1) = R(k,1) + dt * D * dyf / 2.0 * u(i,j-1,0)
			endif

		enddo
	enddo
	
!	Call LAPACK routine to solve backsubstitution
	call SPOTRS( 'U', size(a,1), 1, A, size(a,1), R(:,1), size(a,1), INFO )
	if (INFO .ne. 0) write(*,*) 'ERROR CALCULATING SOLUTION FROM CHOLESKY DECOMPOSITION.... FAILED'

	! For debug - write matrix w (solution to Ax=B)
	!call writevector(w,filename2)
	!call writevector(R(:,1),filename2)

	! Recreate u from w
	do i=1,bigI
		do j=1,bigJ
			n = i + (j-1)*bigI
			u(i,j,0) = R(n,1) !w(n)
		enddo
	enddo
	
END SUBROUTINE CRANK_NICHOLSON_2D_SOLVER

! Subroutine CN_FTCS_FULL_SOLVER is a Crank-Nicolson - FTCS solver for the advection-diffusion equation
! 
! ** Input parameters **
! x : input x vector (evenly spaced)
! y : input y vector (evenly spaced)
! bigI : input mesh I value
! bigJ : input mesh J value
! dt : timestep
! D : diffusion coefficient
! fg1,fg2,fg3,fg4 : dirichlet boundary conditions
! f3a,f3b : advection boundary functions
! A : 'U' cholesky decomposed matrix of LHS (n+1 timestep)
! u : output solution mesh matrix
SUBROUTINE CN_FTCS_FULL_SOLVER(x,y,bigI,bigJ,dt,D,fg1,fg2,fg3,fg4,f3a,f3b,A,u)
	use LinAl
	IMPLICIT NONE
	integer :: i,j,k,m,n,bigI,bigJ,INFO
	real :: dy,dx,dxf,dyf,dt,D!,C
	real, dimension(0:bigI+1) :: x
	real, dimension(0:bigJ+1) :: y
	real, dimension(1:bigI*bigJ, 1:bigI*bigJ) :: A
	real, dimension(bigI*bigJ,1) :: R ! B
	real, dimension(0:bigI+1,0:bigJ+1,0:2) :: u
	real, dimension(bigI*bigJ) :: w
	character (len=100) :: filename = 'hw6debugA.dat'   ! Output data filename
	character (len=100) :: filename2 = 'hw6debugB.dat'   ! Output data filename
	real :: fg1,fg2,fg3,fg4,f3a,f3b
	external fg1,fg2,fg3,fg4,f3a,f3b
	
	R = 0.0
	
	! Calculate dy, dx (assume x and y vectors are evenly spaced)
	dy = y(1) - y(0)
	dx = x(1) - x(0)
	dyf = 1. / (dy ** 2)
	dxf = 1. / (dx ** 2)
	
	do i=1,bigI
		do j=1,bigJ
			k = i + (j-1) * bigI
			
			R(k,1) = u(i,j,0) - u(i,j,0) * dt * D * (dxf + dyf) !  CN Contribution

			! Check for indexing to boundary conditions
			if (i .eq. bigI) then ! i=bigI
				R(k,1) = R(k,1) + dt * D * dxf * fg4(y(j))   !  CN Contribution
				R(k,1) = R(k,1) - dt * f3a(x(i),y(j)) * dxf * fg4(y(j))   !  FTCS Contribution
			else  ! i < bigI
				R(k,1) = R(k,1) + dt * D * dxf / 2.0 * u(i+1,j,0)   !  CN Contribution
				R(k,1) = R(k,1) - dt * f3a(x(i),y(j)) * dxf / 2.0 * u(i+1,j,0)   ! FTCS contribution
			endif
			
			if (i .eq. 1) then ! i=1
				R(k,1) = R(k,1) + dt * D * dxf  * fg3(y(j))   !  CN Contribution
				R(k,1) = R(k,1) + dt * f3a(x(i),y(j)) * dxf  * fg3(y(j))   !  FTCS Contribution
			else ! i > 1
				R(k,1) = R(k,1) + dt * D * dxf / 2.0 * u(i-1,j,0)   !  CN Contribution
				R(k,1) = R(k,1) + dt * f3a(x(i),y(j)) * dxf / 2.0 * u(i-1,j,0)   ! FTCS contribution
			endif
			
			
			if (j .eq. bigJ) then  ! j=bigJ
				R(k,1) = R(k,1) + dt * D * dyf  * fg2(x(i))   !  CN Contribution
				R(k,1) = R(k,1) - dt * f3b(x(i),y(j)) * dyf  * fg2(x(i))   !  FTCS Contribution
			else  ! j < bigJ
				R(k,1) = R(k,1) + dt * D * dyf / 2.0 * u(i,j+1,0)   !  CN Contribution
				R(k,1) = R(k,1) - dt * f3b(x(i),y(j)) * dyf / 2.0 * u(i,j+1,0)   !  FTCS Contribution
			endif
			
			if (j .eq. 1) then ! j=1
				R(k,1) = R(k,1) + dt * D * dyf  * fg1(x(i))   !  CN Contribution
				R(k,1) = R(k,1) + dt * f3b(x(i),y(j)) * dyf  * fg1(x(i))   !  FTCS Contribution
			else ! j > 1
				R(k,1) = R(k,1) + dt * D * dyf / 2.0 * u(i,j-1,0)   !  CN Contribution
				R(k,1) = R(k,1) + dt * f3b(x(i),y(j)) * dyf / 2.0 * u(i,j-1,0)   !  FTCS Contribution
			endif

		enddo
	enddo

!	Call LAPACK routine to solve backsubstitution
	call SPOTRS( 'U', size(a,1), 1, A, size(a,1), R(:,1), size(a,1), INFO )
	if (INFO .ne. 0) write(*,*) 'ERROR CALCULATING SOLUTION FROM CHOLESKY DECOMPOSITION.... FAILED'

	! For debug - write matrix w (solution to Ax=B)
	!call writevector(w,filename2)
	!call writevector(R(:,1),filename2)

	! Recreate u from w
	do i=1,bigI
		do j=1,bigJ
			n = i + (j-1)*bigI
			u(i,j,0) = R(n,1) !w(n)
		enddo
	enddo

END SUBROUTINE CN_FTCS_FULL_SOLVER


! Subroutine  FTCS_2D_SOLVER is a FTCS solver for the 2D diffusion equation
! 
! ** Input parameters **
! x : input x vector (evenly spaced)
! y : input y vector (evenly spaced)
! bigI : input mesh I value
! bigJ : input mesh J value
! dt : timestep
! D : diffusion coefficient
! fg1,fg2,fg3,fg4 : dirichlet boundary conditions
! u : output solution mesh matrix
SUBROUTINE FTCS_2D_SOLVER(x,y,bigI,bigJ,dt,D,fg1,fg2,fg3,fg4,u)
	use LinAl
	IMPLICIT NONE
	integer :: i,j,k,m,n,bigI,bigJ,INFO
	real :: dy,dx,dxf,dyf,dt,D
	real, dimension(0:bigI+1) :: x
	real, dimension(0:bigJ+1) :: y
	real, dimension(1:bigI*bigJ, 1:bigI*bigJ) :: A
	real, dimension(bigI*bigJ,1) :: R ! B
	real, dimension(0:bigI+1,0:bigJ+1,0:2) :: u
	real, dimension(bigI*bigJ) :: w
	character (len=100) :: filename = 'hw6debugA.dat'   ! Output data filename
	character (len=100) :: filename2 = 'hw6debugB.dat'   ! Output data filename
	real :: fg1,fg2,fg3,fg4,fsolve
	external fg1,fg2,fg3,fg4,fsolve
	
	! Calculate dy, dx (assume x and y vectors are evenly spaced)
	dy = y(1) - y(0)
	dx = x(1) - x(0)
	dyf = 1. / (dy ** 2)
	dxf = 1. / (dx ** 2)
	
	do i=1,bigI
		do j=1,bigJ
			k = i + (j-1) * bigI
		
			u(i,j,1) = u(i,j,0) * ( 1.0 - 2 * dt * D * (dxf + dyf) ) 

		! Check for indexing to boundary conditions
			if (i .lt. bigI) then
				if (i .gt. 1) then
					u(i,j,1) = u(i,j,1) + dt * D * dxf * u(i+1,j,0)
					u(i,j,1) = u(i,j,1) + dt * D * dxf * u(i-1,j,0)
				else  ! if i=1
					u(i,j,1) = u(i,j,1) + dt * D * dxf * fg3(y(j)) + dt * D * dxf * u(i+1,j,0)
				endif
			else  ! if i=bigI
				u(i,j,1) = u(i,j,1) + dt * D * dxf * fg4(y(j)) + dt * D * dxf * u(i-1,j,0)
			endif
			
			if (j .lt. bigJ) then
				if (j .gt. 1) then
					u(i,j,1) = u(i,j,1) + dt * D * dyf * u(i,j-1,0)
					u(i,j,1) = u(i,j,1) + dt * D * dyf * u(i,j+1,0)
				else ! if j=1
					u(i,j,1) = u(i,j,1) + dt * D * dyf * fg1(x(i)) + dt * D * dyf * u(i,j+1,0)
				endif
			else  ! if j=bigJ
				u(i,j,1) = u(i,j,1) + dt * D * dyf * fg2(x(i)) + dt * D * dyf * u(i,j-1,0)
			endif
			
		enddo
	enddo
	
	u(:,:,0) = u(:,:,1)
	u(:,:,1) = 0.0
	
	! Populate new u with boundaries from dirichlet conditions
	do i=0,bigI+1
		u(i,0,0) = fg1(x(i))
		u(i,bigJ+1,0) = fg2(x(i))
	enddo
	
	do j=0,bigJ+1
		u(0,j,0) = fg3(y(j))
		u(bigI+1,j,0) = fg4(y(j))
	enddo
	
END SUBROUTINE FTCS_2D_SOLVER

SUBROUTINE FTCS_1D_SOLVER(x,y,bigI,bigJ,dt,D,fg1,fg2,fg3,fg4,u)
	use LinAl
	IMPLICIT NONE
	integer :: i,j,k,m,n,bigI,bigJ,INFO
	real :: dy,dx,dxf,dyf,dt,D
	real, dimension(0:bigI+1) :: x
	real, dimension(0:bigJ+1) :: y
	real, dimension(1:bigI*bigJ, 1:bigI*bigJ) :: A
	real, dimension(bigI*bigJ,1) :: R ! B
	real, dimension(0:bigI+1,0:bigJ+1,0:2) :: u
	real, dimension(bigI*bigJ) :: w
	character (len=100) :: filename = 'hw6debugA.dat'   ! Output data filename
	character (len=100) :: filename2 = 'hw6debugB.dat'   ! Output data filename
	real :: fg1,fg2,fg3,fg4,fsolve
	external fg1,fg2,fg3,fg4,fsolve
	
	! Calculate dy, dx (assume x and y vectors are evenly spaced)
	dy = y(1) - y(0)
	dx = x(1) - x(0)
	dyf = 0.0 ! 1. / (dy ** 2)
	dxf = 1. / (dx ** 2)
	
	do i=1,bigI
		do j=0,0!bigJ
			k = i + (j-1) * bigI
		
			u(i,j,1) = u(i,j,0) * ( 1.0 - 2 * dt * D * (dxf + dyf) ) 

		! Check for indexing to boundary conditions
			if (i .lt. bigI) then
				if (i .gt. 1) then
				 	u(i,j,1) = u(i,j,1) + dt * D * dxf * u(i+1,j,0)
				 	u(i,j,1) = u(i,j,1) + dt * D * dxf * u(i-1,j,0)
				else  ! if i=1
					u(i,j,1) = u(i,j,1) + dt * D * dxf * 1.0 + dt * D * dxf * u(i+1,j,0)
				endif
			else  ! if i=bigI
				u(i,j,1) = u(i,j,1) + dt * D * dxf * 0.0 + dt * D * dxf * u(i-1,j,0)				
			endif
		enddo
	enddo
	
	u(:,:,0) = u(:,:,1)
	u(:,:,1) = 0.0
	
END SUBROUTINE FTCS_1D_SOLVER

! Final project : Question 3 initial condition
! function u0
function f3(x,y)
	implicit none
	real :: f3,x,y
	f3 = y
end function f3

! Final project : Question 2 initial condition
! function u0
function f2(x,y)
	implicit none
	real :: f2,x,y
	f2 = y + sin(3 * PI * x) * sin(2 * PI * y)
end function f2

! function f
function f1(x,y)
	implicit none
	real :: f1,x,y
	f1 = exp(-(x-0.5)**2/0.01) * exp(-(y-0.5)**2/0.01)!0 ! sin(PI * x) * sin(PI * y)
end function f1

! Final project : Question 3
! function a
function fa(x,y)
	implicit none
	real :: fa,x,y
	fa = PI * sin(2 * PI * x) * cos(PI * y)
end function fa

! Final project : Question 3
! function b
function fb(x,y)
	implicit none
	real :: fb,x,y
	fb = -2 * PI * cos(2 * PI * x) * sin(PI * y)
end function fb

! boundary function g1 (bottom)
!  u(x,0)
function g1(x)
	implicit none
	real :: g1,x
	g1 = 0
end function g1
! boundary function g2 (top)
!  u(x,1)
function g2(x)
	implicit none
	real :: g2,x
	g2 = 1 !sin(PI * x)
end function g2
! boundary function g3 (left)
!   u(0,y)
function g3(x)
	implicit none
	real :: g3,x
	g3 = x
end function g3
! boundary function g4 (right)
!   u(1,y)
function g4(x)
	implicit none
	real :: g4,x
	g4 = x
end function g4

subroutine printresult(x,y,u,time,nfile)
implicit none
! This routine prints out the solution at a particular timestep. 
   real, dimension(0:,0:) :: u
   real, dimension(:) :: x,y
!   real, dimension(0:nmesh+1) :: utheor
   real :: time!,bmt,mpi
   integer :: nmesh,nfile,i,m
   integer, parameter :: mmax = 500
  ! real, parameter :: pi = dacos(-1.d0)
   character hundred,ten,unit
   
   hundred = char(nfile/100+48)
   ten=char(mod(nfile,100)/10+48)
   unit=char(mod(nfile,10) + 48)
      
   open(10,file='u.'//hundred//ten//unit,status='unknown')
   open(13,file='time.'//hundred//ten//unit,status='unknown')
   
! This calculates the theoretical solution for the 
! initial step function problem with specific values
! For ua = 1, ub = 0, xa = 0 and xb = 1.

 !  do i=0,nmesh+1
 !     utheor(i) = 1.d0 - x(i)    
 !  enddo

 !  do m=1,mmax
 !     mpi = dble(m)*pi
 !     bmt = (-2.d0/mpi) * dcos(mpi/2.d0) * dexp(-mpi**2*kappa*time)
 !     do i=0,nmesh+1
 !        utheor(i) = utheor(i) + bmt*dsin(mpi*x(i))
 !     enddo
 !  enddo

! Writes out the solution and the theoretical solution. 
 !  do i=0,nmesh+1
 !     write(10,*) x(i),u(i),utheor(i)
 !  enddo
	do i=1,size(x)-1
		do j=1,size(y)-1
			! 2D - Splot
			write(10,*) x(i),y(j),u(i,j)
			! 1D - plot
			!write(10,*) x(i),u(i,0)
		enddo
		write(10,*) ''
	enddo
	
   
   write(13,*) time
   
   close(10)
   close(13)
   
 end subroutine printresult

! Poisson solver using dirichlet conditions
! Input parameters
! x : input x vector (evenly spaced)
! y : input y vector (evenly spaced)
! bigI : input I value
! bigJ : input J value
! f : function to solve
! g1,g2,g3,g4 : boundary functions
! u : output solution mesh matrix
!
SUBROUTINE POISSON_SOLVE_DIRICHLET(x,y,bigI,bigJ,fsolve,fg1,fg2,fg3,fg4,u)
	use LinAl
	IMPLICIT NONE
	integer :: i,j,k,m,n,bigI,bigJ,INFO
	real :: dy,dx,dxf,dyf
	real, dimension(0:bigI+1) :: x
	real, dimension(0:bigJ+1) :: y
	real, dimension(1:bigI*bigJ, 1:bigI*bigJ) :: A
	real, dimension(bigI*bigJ,1) :: B
	real, dimension(0:bigI+1,0:bigJ+1) :: u
	real, dimension(bigI*bigJ) :: w
	character (len=100) :: filename = 'hw6debugA.dat'   ! Output data filename
	character (len=100) :: filename2 = 'hw6debugB.dat'   ! Output data filename
	real :: fg1,fg2,fg3,fg4,fsolve
	external fg1,fg2,fg3,fg4,fsolve
	
	A = 0
	B = 0
	
	! Calculate dy, dx (assume x and y vectors are evenly spaced)
	dy = y(1) - y(0)
	dx = x(1) - x(0)
	dyf = 1. / (dy ** 2)
	dxf = 1. / (dx ** 2)
	
	! Create A and RHS vector B
	do j=1,bigJ 
		do i=1,bigI  
			! calculate n index
			n = i + (j-1) * bigI
			
			A(n,n) = -2.0 * (dxf + dyf)
			B(n,1) = fsolve(x(i),y(j))
			
				! Check for indexing to boundary conditions
				if (i .lt. bigI) then
					A(n,n+1) = dxf
				else 
					B(n,1) = B(n,1) - dxf * fg4(y(j))
				endif
				
				if (i .gt. 1) then
					A(n,n-1) = dxf
				else
					B(n,1) = B(n,1) - dxf * fg3(y(j))
				endif
				if (j .lt. bigJ) then
					A(n,n+bigI) = dyf  
				else 
					B(n,1) = B(n,1) - dyf * fg2(x(i))
				endif
				
				if (j .gt. 1) then
					A(n,n-bigI) = dyf
				else
					B(n,1) = B(n,1) - dyf * fg1(x(i))
				endif
	!		endif

		enddo
	enddo

	!  For Debug: Output matrix A and B
	!call writematrix(A,filename)
	!call writematrix(B,filename2)
	
	! Solve Ax = B (cholesky)
    A = -A
	
!	Call LAPACK routine to solve cholesky decomposition
	call SPOTRF( 'U', size(a,1), A, size(a,1), INFO )
	if (INFO .ne. 0) write(*,*) 'ERROR CALCULATING CHOLESKY DECOMPOSITION.... FAILED',INFO
 
	w = -B(:,1)
	
!3. Solve RT y = AT b (forward solve), and Solve Rx = y (backward solve) .
!	Call LAPACK routine to solve backsubstitution
	call SPOTRS( 'U', size(a,1), 1, A, size(a,1), w, size(a,1), INFO )
	if (INFO .ne. 0) write(*,*) 'ERROR CALCULATING SOLUTION FROM CHOLESKY DECOMPOSITION.... FAILED'

	! For debug - write matrix w (solution to Ax=B)
	!call writevector(w,filename2)

	! Recreate u from w
	do i=1,bigI
		do j=1,bigJ
			n = i + (j-1)*bigI
			u(i,j) = w(n)
		enddo
	enddo
	
	do i=0,bigI+1
		u(i,0) = fg1(x(i))
		u(i,bigJ+1) = fg2(x(i))
	enddo
	
	do j=0,bigJ+1
		u(0,j) = fg3(y(j))
		u(bigI+1,j) = fg4(y(j))
	enddo
	
END SUBROUTINE POISSON_SOLVE_DIRICHLET

! Write .dat file for gnuplot splot
SUBROUTINE WriteSPlotDat(x,y,a,filename)
	IMPLICIT NONE
	INTEGER :: n,i,j
	CHARACTER(LEN=*) :: filename
	real, dimension(:) :: x,y
	real, dimension(:,:) :: a
	n=18
	open(n,file=filename)

	do i=1,size(x)
		do j=1,size(y)
			write(n,*) x(i),y(j),a(i,j)
		enddo
	enddo
	write(*,*) 'Printing output splot data points to: ',filename
	close(n)
END SUBROUTINE WriteSPlotDat

! deriv is a function deriv(t,y,dy) that takes t and y as inputs and returns dy
!  y and dy are arrays of size (n) where n is the order of the ODE (SET ABOVE AS A PARAMETER)
subroutine deriv(t,y,dy)
	implicit none
	integer :: n
	real :: t!,y
	real, dimension(:) :: y
	real, dimension(:), allocatable :: dy
	
	n = size(y)

!	if (allocated(dy)) then 
!		deallocate(dy)
!	endif
	allocate(dy(n))
	
!    f'' = -f' -f, with initial conditions f(0) = 2, f'(0) = 1. Test different timesteps.
!  solve ordinary differential equation y''(t)=-y'(t)-y(t), y(0)=2, y'(0)=1
!   eigenvalue {{1,0},{-1,-1}}
!   lambda1 = -1,   lambda2 = 1    v1 = (0,1)    v2 = (-2,1)
!     y1 = f,  y2 = f',          y1' = y2   ,   y2' = -y2 - y1
	dy(1) = y(2)
	dy(2) = -y(2) - y(1)

end subroutine deriv

! Set initial conditions for above function
SUBROUTINE initialconditions()
	 !  Initial Condition for function
	 allocate(y0(order))
	y0(1) = 2.
	y0(2) = 1.
END SUBROUTINE initialconditions

!    Create a driver program that integrates the ODE above using one method, and then the other. The program must take in as
!     arguments the initial time, the initial value of f, the final time, and the timestep. The user must know how to modify 
!     the routine that calculates G(f,t).
! Input parameters
!   deriv : derivative function
!   orderf : order of derivative function
!   t0 : initial time
!   tf : final time
!   dt : timestep
!   filename : output data filename for coordinates
!   algorithm : implementation step routine, either 'EulerModified' or 'EulerExplicit'
!
SUBROUTINE ODEDRIVER(derivf,t0,tf,dt,filename,algorithm)
	IMPLICIT NONE
	integer :: i,j,k,m,n
	real :: t0,dt,tf
	external derivf
	type(state) :: curstate
	character (len=*) :: filename,algorithm

	write(*,*) 'Input Values :: timestep: ',dt,'Initial Time: ',t0,'Final Time: ',tf
	if (algorithm .eq. 'EulerExplicit') then
		write(*,*) 'Using Eulers Explicit Algorithm'
	else if (algorithm .eq. 'EulerModified') then
		write(*,*) 'Using Eulers Modified Algorithm'
	else if (algorithm .eq. 'AdamsBashforth2') then
		write(*,*) 'Using Adams-Bashforth 2nd order Algorithm'
	else
		write(*,*) 'Invalid Input Algorithm'
	endif
	
	!  Initial Condition for function
	call initialconditions()
	
	! calculate number of timesteps
	n = tf / dt + 2

	! Initializations
	allocate(curstate%y(order,n))
	do i=1,n
		curstate%y(i,1) = y0(i) 
	enddo
	curstate%t = t0
	curstate%step = 1
	allocate(curstate%x(n))
	curstate%x(1) = t0
	
    do while (curstate%t .lt. tf)
	
		if (algorithm .eq. 'EulerExplicit') then
			call EULEREXPLICIT(curstate,dt,derivf)
		else if ((algorithm .eq. 'EulerModified') .or. ((algorithm .eq. 'AdamsBashforth2') .and. (curstate%step .lt. 3))) then
			call EULERMODIFIED(curstate,dt,derivf)
		else if (algorithm .eq. 'AdamsBashforth2') then
			call ADAMSBASHFORTH(curstate,dt,derivf,2)
		endif

	enddo
	
	! Print output coordinates to file
	open(12,file=filename) 
	do i=1,size(curstate%x)
		write(12,*) curstate%x(i),curstate%y(:,i)
	enddo
	write(*,*) 'Printing output fit coordinates to: ',filename
	close(12)	
	
	deallocate(curstate%y)
	deallocate(curstate%x)
	deallocate(y0)
	
END SUBROUTINE ODEDRIVER


!  Adams-Bashforth 2nd order scheme. 
! the scheme must be able to deal with any number of coupled ODEs. 
!   As in the previous homework: create a timestepper routine, and a driver routine. Make sure the timestepper routine
!   takes in as argument the RHS function/routine that returns the derivatives. Bonus point for elegant use of derived data types.
SUBROUTINE ADAMSBASHFORTH(curstate,dt,derivf,aorder)
	IMPLICIT NONE
	INTEGER :: i,j,k,m,nsize,msize,prevstep,aorder
	real :: f0,t0,t,dt,tf,y0, deriv,x
	real, dimension(:), allocatable :: dy,fn1,fn2
	type(state) :: curstate
	interface
		subroutine derivf(t,y,dy)
			integer :: n
			real :: t
			real, dimension(:) :: y
			real, dimension(:), allocatable :: dy
		end subroutine derivf
	end interface

	prevstep = curstate%step
	curstate%step = curstate%step + 1
	curstate%t = curstate%t + dt
	
	! Calculate f_n-1
	call derivf(curstate%x(prevstep),curstate%y(:,prevstep),fn1)
	! Calculate f_n-2
	call derivf(curstate%x(prevstep-1),curstate%y(:,prevstep-1),fn2)
	! Calculate y(n)
	if (aorder .eq. 2) then	
		curstate%y(:,curstate%step)  = curstate%y(:,prevstep)  + 0.5 * dt * (3.0 * fn1 - fn2)
	endif

	curstate%x(curstate%step) = curstate%t
	
END SUBROUTINE ADAMSBASHFORTH

!    Create a routine that advances the solution to an ODE of the kind df/dt = G(f,t) for one timestep using Euler's modified 
!     midpoint method. The input and outputs to the routine must be a derived data type that contains (1) the time and (2)
!     the value of the function f at that time. The name of the function that contains G(f,t) must be passed as an argument
!     to the routine.
SUBROUTINE EULERMODIFIED(curstate,dt,derivf)
	IMPLICIT NONE
	INTEGER :: i,j,k,m,nsize,msize,prevstep
	real :: f0,t0,t,dt,tf,y0, deriv,x!,k1,k2
	real, dimension(:), allocatable :: dy,k1,k2
	type(state) :: curstate
	interface
		subroutine derivf(t,y,dy)
			integer :: n
			real :: t
			real, dimension(:) :: y
			real, dimension(:), allocatable :: dy
		end subroutine derivf
	end interface

	prevstep = curstate%step
	curstate%step = curstate%step + 1
	curstate%t = curstate%t + dt
	
	allocate(k1(order))
	allocate(k2(order))
	
	nsize = size(curstate%y,1)

	call derivf(curstate%x(prevstep),curstate%y(:,prevstep),dy)
	k1 = dt * dy
	deallocate(dy)
	call derivf(curstate%x(curstate%step),curstate%y(:,prevstep)+k1,dy)
	k2 = dt * dy
	deallocate(dy)
	
	curstate%y(:,curstate%step)  = curstate%y(:,prevstep)  + 0.5 * k1 + 0.5 * k2

	curstate%x(curstate%step) = curstate%t
	
END SUBROUTINE EULERMODIFIED


!    Create a routine that advances the solution to an ODE of the kind df/dt = G(f,t) for one timestep using Euler's method. 
!     The input and outputs to the routine must be a derived data type that contains (1) the time and (2) the value of the
!     function f at that time. The name of the function that contains G(f,t) must be passed as an argument to the routine.
!
!  Approximates y(t) in y'(t) = f(y, t) with y(a) = y0 and t = t0..tf and step size dt.
SUBROUTINE EULEREXPLICIT(curstate,dt,derivf)
	IMPLICIT NONE
	INTEGER :: i,j,k,n,m,nsize,msize,prevstep
	real :: f0,t0,t,dt,tf,y0,x
	real, dimension(:), allocatable :: dy
	type(state) :: curstate
	interface
		subroutine derivf(t,y,dy)
			integer :: n
			real :: t
			real, dimension(:) :: y
			real, dimension(:), allocatable :: dy
		end subroutine derivf
	end interface

	prevstep = curstate%step
	curstate%step = curstate%step + 1
	curstate%t = curstate%t + dt

	nsize = size(curstate%y,1)

	call derivf(curstate%x(prevstep),curstate%y(:,prevstep),dy)

	do i=1,nsize
		curstate%y(i,curstate%step)  =  curstate%y(i,prevstep)  + dt * dy(i)
		curstate%x(curstate%step) = curstate%t
	enddo

END SUBROUTINE EULEREXPLICIT

end module DiffEq
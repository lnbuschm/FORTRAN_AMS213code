!  ********** HOMEWORK #6 --------- AMS 213 *************** !
!  Luke Buschmann
program hw6
	use LinAl
	use DiffEq
	implicit none
	integer :: i,j,k
	integer, parameter :: bigI = 3  ! Input parameter bigI
	integer, parameter :: bigJ = 3  ! Input parameter bigJ
	real, parameter :: dx = 1. / (bigI+1.)
	real, parameter :: dy = 1. / (bigJ+1.) 
	character (len=100) :: filename = 'hw6mesh.dat'   ! Output data filename
	real, dimension(0:bigI+1) :: x
	real, dimension(0:bigJ+1) :: y
	real, dimension(0:bigI+1,0:bigJ+1) :: u
	real, dimension(1:bigI*bigJ, 1:bigI*bigJ) :: A,B
	real :: deriv2,dt
	

	write(*,*) '******** HOMEWORK #6  --  AMS213 ******** '
	
! Create a Poisson-Solver routine that solves Laplacian u = f(x,y) for u, in a rectangular domain,
!  with boundary conditions u(x, y_bottom) = g1(x), u(x,y_top) = u(g2(x), u(x_left,y) = g3(y) and u(x_right,y) = g4(y).
!  	The routine must take as argument the mesh x(0:I+1) and y(0:J+1), the dimensions (if you think it's necessary),
!   the functions f, and g1, g2, g3, g4 (declared as external) and finally, return u(x,y) as an array. 
!   Write a driver program that creates the mesh, calls this routine, and returns the solution.
	
	! Create X and Y vectors
	! Set size of mesh  x=[0,1], y=[0,1]
	x(0) = 0
	x(bigI+1) = 1
	y(0) = 0
	y(bigJ+1) = 1

	write(*,*) 'dx= ',dx,'    dy= ',dy

	! Populate X and Y vectors
	do i=1,bigI
		x(i) = i * dx
	enddo
	
	
	do j=1,bigJ
		y(j) = j * dy
	enddo
	
	write(*,*) 'x vector',x
	write(*,*) 'y vector',y
	
	! Call solver
	write(*,*) 'Calling DIRICHLET POISSON SOLVER...'
	call POISSON_SOLVE_DIRICHLET(x,y,bigI,bigJ,f1,g1,g2,g3,g4,u)
		
	! Write solution mesh to file
	call WriteMatrix(u,filename)
	

end program hw6
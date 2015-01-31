!  ********** FINAL PROJECT --------- AMS 213 *************** !
!  Luke Buschmann
program projectfinal
	use LinAl
	use DiffEq
	implicit none
	integer :: i,j,k,niter,it,nfile,inputOption
	integer, parameter :: bigI = 40 !16  ! Input parameter bigI
	integer, parameter :: bigJ = 40 !16  ! Input parameter bigJ
	integer, parameter :: nOutputFiles = 10  ! # of output files
	real, parameter :: dx = 1. / (bigI+1.)
	real, parameter :: dy = 1. / (bigJ+1.) 
	character (len=100) :: filename = 'hw6mesh.dat'   ! Output u matrix filename
	character (len=100) :: filename2 = 'heat.dat'   ! Output splot data filename
	real, dimension(0:bigI+1) :: x
	real, dimension(0:bigJ+1) :: y
	real, dimension(0:bigI+1,0:bigJ+1,0:2) :: u
	real, dimension(1:bigI*bigJ, 1:bigI*bigJ) :: A,B
	real :: deriv2,dt,t0,tf,t,D,ditprint,diter

	write(*,*) '******** FINAL PROJECT  --  AMS213 ******** '
	
	! Create X and Y vectors
	! Set size of mesh  x=[0,1], y=[0,1]
	x(0) = 0
	x(bigI+1) = 1
	y(0) = 0
	y(bigJ+1) = 1

	! Populate X and Y vectors
	do i=1,bigI
		x(i) = i * dx
	enddo
	
	do j=1,bigJ
		y(j) = j * dy
	enddo

	write(*,*) 'PDE Solvers:  Enter 1 for QUESTION 2: Crank Nicholson, 2 for QUESTION 2: FTCS, '
	write(*,*) '     enter 3 for QUESTION 3: FTCS-advection CN-diffusion'
	read(*,*) inputOption
	
	! Diffusion coefficient D:  Higher value = faster diffusion
	D = 1.0 ! 0.1 !1.0 !0.01  
	t0 = 0.0
	tf = 1.0
	dt = 0.00001
	niter = int((tf-t0)/dt)
	write(*,*) 'Total number of steps: niter =', niter
	 
! *** Enter printing options
!	write(*,*) 'Elapsed time between two outputs?'
!	read(*,*) dtprint
!	write(*,*) 'Initial file number?'
!	read(*,*) nfile 
	nfile = 1 ! Initial file number
	ditprint = niter / nOutputFiles !* dt ! (tf - t0) / 8.0  !  output mesh for ten intervals between t0 and tf
	diter = 1 ! ditprint 

	
	
	! Initialize u at time 0
	u = 0.0
	do i=0,bigI+1
		u(i,0,0) = g1(x(i))
		u(i,bigJ+1,0) = g2(x(i))
	enddo
	
	do j=0,bigJ+1
		u(0,j,0) = g3(y(j))
		u(bigI+1,j,0) = g4(y(j))
	enddo
	do i=1,bigI
		do j=1,bigJ
			! Create u0 for Question 2 (if specified)
			if (inputOption .lt. 3) then
				u(i,j,0) = f2(x(i),y(j))
			else
			! Create u0 for Question 3 (if specified)
				u(i,j,0) = f3(x(i),y(j))
			endif
		enddo
	enddo
	
	write(*,*) 'Outputting ',nOutputFiles,' splot files over time interval t0=',t0,', tf=',tf
	write(*,*) 'Diffusion coefficient D = ',D,' timestep dt = ',dt
	write(*,*) 'dx = ',dx,'    dy = ',dy
		
	t = t0
	
	! For 1D FTCS - initialize u as step function
!	u = 0.0
!	  u(0,0,0) = 1.0
!	  do i=1,bigI
!		 if(x(i).lt.(1.0)/2.0) then
!			u(i,0,0) = 1.0
!		 else
!			u(i,0,0) = 0.0
!		 endif
!	  enddo
!	  u(bigI+1,0,0) = 0.0
	  
	  ! Print initial u
	call printresult(x,y,u(:,:,0),t,nfile)  
	write(*,*) 'wrote file at t=',t !,'u(8,0,0)=',u(50,0,0)
!	write(*,*) 'x vector',x
!	write(*,*) 'y vector',y
	write(*,*) 'Calling solver INIT...'
	! Creates A and performs a cholesky decomposition
	if ((inputOption .eq. 1) .or. (inputOption .eq. 3)) then
		call CRANK_NICHOLSON_2D_INIT(x,y,bigI,bigJ,dt,D,g1,g2,g3,g4,A)
	endif
	
	write(*,*) 'Calling equation SOLVER...'
		
!	do while (t .lt. tf)
! *** Time-step here
	do it = 1,niter
		t = t + dt
		! Call solver
		if (inputOption .eq. 1) then 
			call CRANK_NICHOLSON_2D_SOLVER(x,y,bigI,bigJ,dt,D,g1,g2,g3,g4,A,u)
		elseif (inputOption .eq. 2) then
			call FTCS_2D_SOLVER(x,y,bigI,bigJ,dt,D,g1,g2,g3,g4,u)
		elseif (inputOption .eq. 3) then
			call CN_FTCS_FULL_SOLVER(x,y,bigI,bigJ,dt,D,g1,g2,g3,g4,fa,fb,A,u)
		endif
	!	call FTCS_1D_SOLVER(x,y,bigI,bigJ,dt,D,g1,g2,g3,g4,u)
	!	write(*,*) 'dtprint-1.d-16=',dtprint-1.d-16,'ddt=',ddt,'dtprint=',dtprint
	!	write(*,*) 'wrote file at t=',t,'u(8,0,0)=',u(50,0,0)

		! Print if it is time to print.
		if(diter.ge.ditprint) then
			
			nfile = nfile + 1
			call printresult(x,y,u(:,:,0),t,nfile)
			diter = 1
			write(*,*) 'wrote file at t=',t !,'u(8,0,0)=',u(50,0,0)
		else
			diter = diter + 1
		endif
	enddo
	! Write solution mesh to file
	!call WriteMatrix(u,filename)
	
!	call WriteSPlotDat(x,y,u(:,:,1),filename2)
	

end program projectfinal
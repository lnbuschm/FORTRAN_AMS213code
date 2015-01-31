module LinAl
implicit none
integer :: nsize,msize
integer,private :: i,j,k


contains

!  **************************************** MIDTERM PROJECT ************************************ !
!  **************************************************************************************** !

! Question 1 : QR Decomposition
! Returns QR decomposition of a non-square nxm matrix 
! Parameters :: a(input) nxm matrix, q(output), r(output) 
SUBROUTINE QRDECOMP(a,q,r)
	IMPLICIT NONE
	INTEGER :: i,j,k,m,nsize,msize
	real :: scalef, total=0, sum=0, tau, sigma, sk
	real, dimension(:,:) :: a
	real, dimension(:,:), allocatable :: Identity,qt,r,rinv,w,P,q
	nsize = size(a,1)
	msize = size(a,2)
	
	! Create identity matrix
	allocate(Identity(nsize,nsize))  ! Identity(nsize,msize))  
	allocate(w(nsize,1))
	Identity=0
	do i=1,nsize
		Identity(i,i)= 1
	enddo
	if (nsize .le. msize) then
		allocate(qt(nsize,nsize))  ! qt is nxn
		allocate(q(nsize,nsize))  ! qt is nxn
		allocate(r(nsize,msize))  ! R is nxn
	else
		allocate(qt(nsize,msize))  ! qt is nxm
		allocate(q(nsize,msize))  ! qt is nxm
		allocate(r(msize,msize))  ! R is nxn
	end if
	
	! Initialize Q as the identity matrix
	qt=Identity 

	! Since A becomes R, we will copy A into R to use in our algorithm
	R = A	
	w = 0
	
	do k=1,msize  !nsize-1
	
		! Calculate norm : sqrt of sum of squares of column k 
		sum = 0.
		do i=k,nsize
			sum = sum + r(i,k) ** 2
		enddo
		
		! ensure sk - a(k,k) != 0
		sk = sqrt(sum)  
		if (sk.eq.r(k,k)) sk = -sk

		! Create w
		do i=1,nsize
			if (i.lt.k) then
				w(i,1) = 0
			else if (i.eq.k) then
				w(i,1) = (r(i,k)-sk) / sqrt(2 * sk * (sk - r(k,k)))
			else 
				w(i,1) = r(i,k) / sqrt(2 * sk * (sk - r(k,k)))
			endif
		enddo
		
		! Create householder matrix P
		P = Identity - 2 * matmul(w,transpose(w))
		
		! Iterate Pk * A
		R = matmul(P,R)

		! zero out elements that are near-zero due to machine accuracy to create a true upper triangular matrix R
		do i=k+1,nsize
			r(i,k) = 0
		enddo
		
		! Iterate Qk = P * Qk
		qt = matmul(P,qt)
		
	enddo
	
	q = transpose(qt)
	qt = transpose(qt)
	
	! Optional:  Resize matrix Q since it will be padded with 0s after transposing
!	if (nsize .lt. msize) then
!		deallocate(q)
!		allocate(q(nsize,nsize))  ! q is nxn
!		do j=1,nsize
!			do i=1,nsize
!				q(i,j) = qt(i,j)
!			enddo
!		enddo
!	else
!		deallocate(q)
!		allocate(q(nsize,msize))  ! q is nxn
!		do j=1,msize
!			do i=1,nsize
!				q(i,j) = qt(i,j)
!			enddo
!		enddo		
!	end if
		
!	write(*,*) '********  Output Matrix Q : QR Decomposition of A ******** '
!	call OutputMatrix(Q)
!	write(*,*) '********  Output Matrix R : QR Decomposition of A ******** '
!	call OutputMatrix(R)
	
	! Check result by multiplying Q*R and verifying it is equal to A (to numerical precision)
!	write(*,*) '********  Output Matrix Q*R: QR Decomposition of A ******** '
!	A = matmul(Q,R)
!	call OutputMatrix(A)
	
	
END SUBROUTINE QRDECOMP


! Question 2 : QR Backsubstitution
! Solves Rm x =  Qtm b
!  Where b is RHS vector of size mx1 and where Rm and Qtm contain only the first m lines of
!    Qt and R matrices obtained as a result of Question 1
!
! Solve A x = B   using  Q R x = b
! Input: nxn matrix A and nx1 matrix B,  Output: Solution in matrix B
SUBROUTINE QRSOLVE(a,b)
	IMPLICIT NONE
	INTEGER :: i,j,k,m,nsize,msize,rhsVecCount
	real :: scalef, total=0, sum=0, tau, sigma
	real, dimension(:,:) :: a, b
	real, dimension(:,:), allocatable :: q,r,qt,qtb
	nsize = size(a,1)
	msize = size(a,2)
	
	! Multiple RHS vectors is not yet implemented, so set k=1
	rhsVecCount = 1
	k = 1
	
	! Obtain QR decomposition of A
	call qrdecomp(a,q,r)

! Solve QR x = b
	! First calculate  R x = Qt b
	
	! Calculate Qt * b
	qt = transpose(q)
	qtb = matmul(qt,b)

	! back substitution of  R x = Qt b
	do i=nsize,1,-1
	!	do k=1,rhsVecCount
			sum = qtb(i,k)
			do j=i+1,nsize
				sum = sum - r(i,j) * qtb(j,k)
			enddo	
			qtb(i,k) = sum / r(i,i)
	!	enddo
	enddo

	! Return solution
	b = qtb
	
END SUBROUTINE QRSOLVE

! Question 3 : Comparison of QR and Cholesky with Atkinson data
! Read atkinson.dat and produce cubic least square fit using QR method
! Returns parameters.dat that returns parameters of the first and RMS error of the fit
! Returns fitted data in fitted.dat formatted as x,y coordinates
SUBROUTINE QRLEASTSQUARES(order,parameters,filename,filename2,filename3)
	IMPLICIT NONE 
	INTEGER :: NR, ios, nsize,msize,i,j,k,order
	real :: sum,ax,bx,cx,error
	real, dimension(:,:), allocatable :: a,aa,bb,q,r,b,qtb,qt,y
	real, dimension(:), allocatable :: parameters
	INTEGER, PARAMETER :: maxrecs = 10000 
	character (len=*) :: filename,filename2,filename3
	CHARACTER(1) :: junk 

	!Determine total number of lines in file 
	NR = 0 
	open(unit=1,file=filename) 
	do j=1,maxrecs 
		read(1,*,iostat=ios) junk 
		if (ios /= 0) exit 
		if (j == maxrecs) then 
			write(*,*) 'error: maximum number of records exceeded...' 
			write(*,*) 'exiting program now...' 
			stop 
		endif 
		NR = NR + 1 
	enddo 
	rewind(1) 

	write(*,*) 'lines in input file,',filename,' =',nr
	
	! allocate data variables 
	allocate(a(nr,1+order))
	allocate(b(nr,1))
	allocate(aa(1+order,1+order))
	allocate(bb(1+order,1))
	
	! Read coordinates from file
	do i=1,NR
		a(i,1)=1.
		read(1,*) a(i,2),b(i,1)
		do k=3,order+1
			a(i,k) = a(i,2) ** (k-1)
		enddo
	enddo
	close(1)

	! Save B matrix ( y values) to calculate error later
	y = b
	
	nsize = size(a,1)
	msize = size(a,2)

! Obtain QR decomposition of A
	call qrdecomp(a,q,r)
	
	! Calculate Qt * b
	q = transpose(q)

	! resize Q
	qt = q  ! use qt as temp storage
	
	if (nsize .lt. msize) then
		deallocate(q)
		allocate(q(nsize,nsize))  ! q is nxn
		do j=1,nsize
			do i=1,nsize
				q(i,j) = qt(i,j)
			enddo
		enddo
	else
		deallocate(q)
		allocate(q(msize,nsize))  ! q is nxn
		do j=1,nsize
			do i=1,msize
				q(i,j) = qt(i,j)
			enddo
		enddo
		
	end if

	! Solve QR x = b
	! First calculate  R x = Qt b
	qtb = matmul(q,b)
	
	! back substitution of  R x = Qt b
	do i=msize,1,-1
	!	do k=1,rhsVecCount
			sum = qtb(i,1)
			do j=i+1,msize
				sum = sum - r(i,j) * qtb(j,1)
			enddo	
			qtb(i,1) = sum / r(i,i)
	!	enddo
	enddo

	! Return solution
	b = qtb

	write(*,*) '********  Least-square coefficient values (QR) ******** '
	call outputmatrix(b)

	! Print output coordinates and parameters to file
	allocate(parameters(size(a,2)))
	open(2,file=filename2) ! Parameter file
	open(3,file=filename3) ! Coordinate file
	error = 0.
	do i=1,size(a,1)
		sum = 0.
		do k=1,order
			sum = sum + b(k+1,1) * a(i,2) ** k
			
		enddo
		sum = sum + b(1,1)
		write(3,*) a(i,2),sum
		error = error + (sum - y(i,1)) ** 2

	enddo
	write(*,*) 'Printing output fit coordinates to:',filename3
	close(3)	
	
	write(2,*) 'Fit Parameters for Polynomial of Order:,',order
	do i=1,size(b,1)
		write(2,*) b(i,1)
	enddo
	
	error = sqrt(error/size(a,1))
	write(2,*) 'RMS Error:',error
	write(*,*) 'Printing output fit parameters to:',filename2
	close(2)	
	
END SUBROUTINE QRLEASTSQUARES


! Question 4 : Failure of Cholesky Method
! Modify QR fit and Cholesky fit algorithms to include an input parameter to set the order of the fit polynomial
!    Write a program that uses the Cholesky decomposition routines you just wrote to find a cubic 
! Least-Square fit to the data given in the Atkinson textbook. The data is tabulated in this file. 
! Print the results into a file, for the fitted parameters, and the fitting curve (as shown in class).
SUBROUTINE CHOLESKYLEASTSQUARES(order,parameters,filename,filename2,filename3)
	IMPLICIT NONE 
	INTEGER :: NR, ios, nsize,msize,i,j,k,order
	real :: sum,error
	real, dimension(:,:), allocatable :: a,b,aT,r,rT,l,aa,bb,y
	real, dimension(:), allocatable :: parameters
	INTEGER, PARAMETER :: maxrecs = 10000 
	character (len=*) :: filename,filename2,filename3
	CHARACTER(1) :: junk 

	!Determine total number of lines in file 
	NR = 0 
	open(unit=1,file=filename) 
	do j=1,maxrecs 
		read(1,*,iostat=ios) junk 
		if (ios /= 0) exit 
		if (j == maxrecs) then 
			write(*,*) 'error: maximum number of records exceeded...' 
			write(*,*) 'exiting program now...' 
			stop 
		endif 
		NR = NR + 1 
	enddo 
	rewind(1) 

	write(*,*) 'lines in input file,',filename,' =',nr
	
	! allocate data variables 
	allocate(a(nr,1+order))
	allocate(b(nr,1))
	allocate(aa(1+order,1+order))
	allocate(bb(1+order,1))
	
	! Read coordinates from file
	do i=1,NR
		a(i,1)=1.
		read(1,*) a(i,2),b(i,1)
		do k=3,order+1
			a(i,k) = a(i,2) ** (k-1)
		enddo
	enddo
	close(1)

	! Save B matrix ( y values) to calculate error later
	y = b
	
!1.	Compute AT A and AT b
	aT = TRANSPOSE (a)
	aa = matmul(aT,a)
	bb = matmul(aT,b)

!2. Compute the Cholesky-decomposition     AT A = RT R.
	call CHOLESKYD(aa,l)
!	write(*,*) '********  Output Matrix L : Cholesky Decomposition of aT * a ******** '
!	call OutputMatrix(l)
	nsize = size(l,1)
	
!3. Solve RT y = AT b (forward solve)
	! forward substitution L y = b
	do i=1,nsize
			sum = bb(i,1)
			do j=i-1,1,-1
				sum = sum - l(i,j) * bb(j,1)
			enddo
			bb(i,1) = sum / l(i,i)
	enddo
	
!4.	Solve Rx = y (backward solve) .
	! back substitution L^T x = y
	do i=nsize,1,-1
			sum = bb(i,1)
			do j=i+1,nsize
				sum = sum - l(j,i) * bb(j,1)
			enddo	
			bb(i,1) = sum / l(i,i)
	enddo
			
	write(*,*) '********  Least-square coefficient values (Cholesky) ******** '
	call OutputMatrix(bb)

	! Print output coordinates and parameters to file
	allocate(parameters(size(a,2)))
	open(2,file=filename2) ! Parameter file
	open(3,file=filename3)
	error = 0.
	do i=1,size(a,1)
		sum = 0.
		do k=1,order
			sum = sum + bb(k+1,1) * a(i,2) ** k
		enddo
		sum = sum + bb(1,1)
		write(3,*) a(i,2),sum
		error = error + (sum - y(i,1)) ** 2
	enddo
	write(*,*) 'Printing output fit coordinates to:',filename3
	close(3)	
	
	write(2,*) 'Fit Parameters for Polynomial of Order:,',order
	do i=1,size(bb,1)
		write(2,*) bb(i,1)
	enddo
	
	error = sqrt(error/size(a,1))
	write(2,*) 'RMS Error:',error
	write(*,*) 'Printing output fit parameters to:',filename2
	close(2)	
	
	
END SUBROUTINE CHOLESKYLEASTSQUARES






!  **************************************** HOMEWORK 3 ************************************ !
!  **************************************************************************************** !
!    Write a routine that takes a real symmetric matrix and performs a Cholesky decomposition. 
! Write another one that takes the output of the Cholesky decomposition, and returns the solution 
! to the problem AX = B where B is a matrix that contains many RHS vectors. Alternatively, just 
! write one routine that does both the Cholesky decomposition and returns the solution to AX=B. 
! Make sure the code is written in such a way that it is possibly to *not* give a B matrix, if 
! the user just wants to do a Cholesky decomposition.
!
!SUBROUTINE CHOLESKYD(a,l)  ::   Input n,n matrix a  ::  Outputs Cholesky decomposition n,n matrix l
SUBROUTINE CHOLESKYD(a,l)
	IMPLICIT NONE
	INTEGER :: i,j,k,m,nsize,msize
	real :: scalef, total=0, sum=0
	real, dimension(:,:) :: a
	real, dimension(:,:), allocatable :: l
	nsize = size(a,1)
	msize = size(a,2)
	allocate(l(nsize,msize))
	l=0
	do j=1,nsize
		total = 0.
		do k=1,j-1
			total = total + (l(j,k) ** 2)
		enddo
		l(j,j) = sqrt(a(j,j) - total)
	!	write(*,*) 'cholesky loop....  j=',j
		do i=j+1,nsize
			total = 0.
			do k=1,j-1
				total = total + l(i,k) * l(j,k)
			enddo
			l(i,j) = ( a(i,j) - total ) / l(j,j) 
		enddo
	enddo

END SUBROUTINE CHOLESKYD

!SUBROUTINE CHOLESKYD(a,l)  ::   Input nxn matrix a  ::  Outputs b matrix with solution
SUBROUTINE CHOLESKYDSOL(a,b)
	IMPLICIT NONE
	INTEGER :: i,j,k,m,nsize,msize,rhsVecCount
	real :: scalef, total=0, sum=0
	real, dimension(:,:) :: a,b
	real, dimension(:,:), allocatable :: l, sol
	
	nsize = size(a,1)
	msize = size(a,2)
	rhsVecCount = size(b,2)
	allocate(sol(nsize,rhsVecCount))
	sol=0
	write(*,*) ' cholesky decomp'
	! perform cholesky decomposition on input matrix A
	call CHOLESKYD(a,l)
	write(*,*) 'Performed cholesky decomp'
	
	! forward substitution L y = b
	do i=1,nsize
		do k=1,rhsVecCount
			sum = b(i,k)
			do j=i-1,1,-1
				sum = sum-l(i,j)*b(j,k)
			enddo
			b(i,k)=sum / l(i,i)
		enddo
	enddo
	
	! back substitution L^T x = y
	do i=nsize,1,-1
		do k=1,rhsVecCount
			sum = b(i,k)
			do j=i+1,nsize
				sum = sum - l(j,i) * b(j,k)
			enddo
			
			b(i,k) = sum / l(i,i)
		enddo
	enddo
		
END SUBROUTINE CHOLESKYDSOL

!    Write a program that uses the Cholesky decomposition routines you just wrote to find a cubic 
! Least-Square fit to the data given in the Atkinson textbook. The data is tabulated in this file. 
! Print the results into a file, for the fitted parameters, and the fitting curve (as shown in class).
SUBROUTINE CUBICLEASTSQUARES(ax,bx,cx,filename,filename3)
	IMPLICIT NONE 
	INTEGER :: NR, ios, nsize,msize,i,j,k
	real :: sum,ax,bx,cx
	real, dimension(:,:), allocatable :: a,b,aT,r,rT,l,aa,bb
	INTEGER, PARAMETER :: maxrecs = 10000 , order = 3 ! 3 for cubic. 2 for quadratic
	character (len=*) :: filename,filename3
	CHARACTER(1) :: junk 

	!Determine total number of lines in file 
	NR = 0 
	open(unit=1,file=filename) 
	do j=1,maxrecs 
		read(1,*,iostat=ios) junk 
		if (ios /= 0) exit 
		if (j == maxrecs) then 
			write(*,*) 'error: maximum number of records exceeded...' 
			write(*,*) 'exiting program now...' 
			stop 
		endif 
		NR = NR + 1 
	enddo 
	rewind(1) 

	write(*,*) 'lines in input file =',nr
	
	! allocate data variables 
	allocate(a(nr,1+order))
	allocate(b(nr,1))
	allocate(aa(1+order,1+order))
	allocate(bb(1+order,1))
	
	! Read coordinates from file
	do i=1,NR
		a(i,1)=1.
		read(1,*) a(i,2),b(i,1)
		do k=3,order+1
			a(i,k) = a(i,2) ** (k-1)
		enddo
	enddo
	close(1)

!1.	Compute AT A and AT b
	aT = TRANSPOSE (a)
	aa = matmul(aT,a)
	bb = matmul(aT,b)
!	write(*,*) '********  Output Matrix AA :  ******** '
!	call OutputMatrix(aa)
	
!2. Compute the Cholesky-decomposition     AT A = RT R.
	call CHOLESKYD(aa,l)
!	write(*,*) '********  Output Matrix L : Cholesky Decomposition of aT * a ******** '
!	call OutputMatrix(l)
	nsize = size(l,1)
	
!3. Solve RT y = AT b (forward solve)
	! forward substitution L y = b
	do i=1,nsize
			sum = bb(i,1)
			do j=i-1,1,-1
				sum = sum - l(i,j) * bb(j,1)
			enddo
			bb(i,1) = sum / l(i,i)
	enddo
	
!4.	Solve Rx = y (backward solve) .
	! back substitution L^T x = y
	do i=nsize,1,-1
			sum = bb(i,1)
			do j=i+1,nsize
				sum = sum - l(j,i) * bb(j,1)
			enddo	
			bb(i,1) = sum / l(i,i)
	enddo
			
	write(*,*) '********  Matrix BB: Cubic least-square coefficient values (HOMEMADE ROUTINE) ******** '
	call OutputMatrix(bb)

	! Print output coordinates to file
	open(3,file=filename3)
	do i=1,size(a,1)
		write(3,*) a(i,2),bb(4,1)*a(i,2)**3 + bb(3,1)*a(i,2)**2 + bb(2,1)*a(i,2) + bb(1,1)
	enddo
	write(*,*) 'Printing output fit coordinates to:',filename3
	close(3)	
	
END SUBROUTINE CUBICLEASTSQUARES

!    Study the LAPACK routines repository to find the correct Cholesky decomposition and backsusbtitution
!  routines for your purposes. Be careful, there are a whole bunch of "Cholesky" routines, but only one is
!  the right one for general-purpose real symmetric positive-definite matrices.
!    Write a second program that uses these LAPACK routines instead to solve the same Least-Square problem 
! as above. The code should be quite similar to your previous one, aside from the routine calls and a few 
! more variable declarations/initializations. Create a Makefile to compile this code and link it to the LAPACK
!  libraries. Run the code, and print the outputs in new files. Compare the outputs of the two codes using your
!  routine, and the LAPACK routines.

SUBROUTINE LAPACKLEASTSQUARES(ax,bx,cx,filename,filename2)
	IMPLICIT NONE 
	INTEGER :: NR, ios, nsize,msize,i,j,k, INFO, rhsVecCount
	REAL :: sum,ax,bx,cx
	REAL, dimension(:,:), allocatable :: a,b,aT,r,rT,l,aa,bb
	INTEGER, PARAMETER :: maxrecs = 10000 , offset = 2 ! 2 for cubic. 1 for quadratic
	character (len=*) :: filename,filename2
	CHARACTER(1) :: junk, UPLO

	!Determine total number of lines in file 
	NR = 0 
	open(unit=1,file=filename) 
	do j=1,maxrecs 
		read(1,*,iostat=ios) junk 
		if (ios /= 0) exit 
		if (j == maxrecs) then 
			write(*,*) 'error: maximum number of records exceeded...' 
			write(*,*) 'exiting program now...' 
			stop 
		endif 
		NR = NR + 1 
	enddo 
	rewind(1) 
	
	write(*,*) 'Lines in Input File =',NR
	
	! Allocate data variables 
	allocate(a(NR,2+offset))
	allocate(b(NR,1))
	allocate(aa(2+offset,2+offset))
	allocate(bb(2+offset,1))
	
	! Read coordinates from file
	do i=1,NR
		a(i,1)=1.
		read(1,*) a(i,2),b(i,1)
		do k=3,offset+2
			a(i,k) = a(i,2) ** (k-1)
		enddo
	enddo

	close(1)
	
!1.	Compute AT A and AT b
	aT = transpose (a)
	aa = matmul(aT,a)
	bb = matmul(aT,b)

!2. Compute the Cholesky-decomposition     AT A = RT R.
	UPLO = 'U'
	nsize = size(aa,1)
	rhsVecCount = 1 ! size(b,2)

!	Call LAPACK routine to solve cholesky decomposition
	call SPOTRF( UPLO, nsize, aa, nsize, INFO )
	if (INFO .ne. 0) write(*,*) 'ERROR CALCULATING CHOLESKY DECOMPOSITION.... FAILED',INFO
 
!	write(*,*) '********  Output Matrix L : Cholesky Decomposition of aT * a ******** '
!	call OutputMatrix(aa)
	
!3. Solve RT y = AT b (forward solve), and Solve Rx = y (backward solve) .
!	Call LAPACK routine to solve backsubstitution
	call SPOTRS( UPLO, nsize, rhsVecCount, aa, nsize, bb, nsize, INFO )
	if (INFO .ne. 0) write(*,*) 'ERROR CALCULATING SOLUTION FROM CHOLESKY DECOMPOSITION.... FAILED'
 
	write(*,*) '********  Matrix BB: Cubic least-square coefficient values (LAPACK ROUTINE) ******** '
	call OutputMatrix(bb)

	! Print output coordinates to file
	open(2,file=filename2)
	do i=1,size(a,1)
		write(2,*) a(i,2),bb(4,1)*a(i,2)**3 + bb(3,1)*a(i,2)**2 + bb(2,1)*a(i,2) + bb(1,1)
	enddo
	write(*,*) 'Printing output fit coordinates to:',filename2
	close(2)	
	
	
END SUBROUTINE LAPACKLEASTSQUARES



!  **************************************** HOMEWORK 2 ************************************ !

!    Write a routine that performs an LU decomposition, and returns the LU form of the matrix. 
!     The LU form must be printed to the screen by the calling program. Pivoting is not necessary.
!    (Optional, could be useful for you later): Copy the LU decomposition and backsubstitution routines from
!     Numerical Recipes or any other source, and learn how to use it to solve linear algebra problems. 
!     You're free to choose a F77 or an F90 source.
SUBROUTINE LU_Decomposition(a,L,U)
	IMPLICIT NONE
	INTEGER :: i,j,k,m,nsize,msize
	real :: scalef, total=0
	real, dimension(:,:) :: a
	real, dimension(:,:), allocatable :: L,U
	nsize = size(a,1)
	msize = size(a,2)
	allocate(L(nsize,msize))
	allocate(U(nsize,msize))
	L = 0
	U = 0
	do i=1,nsize
		L(i,i) = 1.
	enddo
	! Crout algorithm
	do k=1,nsize
		
		do j=k,nsize
			total=0
			do m=1,k-1
				total = total + L(k,m) *  U(m,j)
			enddo
	!		write (*,*) 'U(k,j) = ',U(k,j),' = ','a(k,j)',a(k,j),'- total',total
			U(k,j) = a(k,j) - total
		enddo
		
		do i=k+1,nsize
			total=0
			do m=1,k-1
				total = total + L(i,m) *  U(m,k)
			enddo
			L(i,k) = (a(i,k) - total) / U(k,k)
		enddo
	enddo


END SUBROUTINE LU_Decomposition

!    Write a routine that returns the inverse of a matrix. The matrix must be read from a file
!    (as in the example given in section), and the result must be printed into a file by the calling program.
SUBROUTINE Gauss_Inverse(a,inv)
	IMPLICIT NONE
	INTEGER :: n=0, pline=1, rhsVecCount,i,j,k
	real :: total,pmax,scalef
	real, dimension(:,:) :: a
	real, dimension(:), allocatable :: matScale
	real, dimension(:,:), allocatable :: inv
	nsize = size(a,1)
	msize = size(a,2)

	allocate(matScale(nsize))
	allocate(inv(nsize,msize))

	!create identity matrix
	do i=1,nsize
		do j=1,msize
			if (i .eq. j) then
				inv(i,j)=1.
			else 
				inv(i,j)=0.
			endif
		enddo
	enddo
	
	! (1) find scale factors for each row
	do i=1,nsize
		do j=1,msize
			if (matScale(i) < a(i,j)) then
				matScale(i) = a(i,j)
			endif
		enddo
	enddo
!	write(*,*) 'scale matrix:'
!	write (*,'(3f8.3)') matScale
	
	
	! (2) scale, pivot, and eliminate
	do j=1,nsize
		! find largest scaled pivot
		pmax=0
		pline=j
		do k=j,nsize
			if (pmax .lt. (a(k,j) / matScale(k))) then
				pmax = a(k,j) / matScale(k)
				pline = k
			endif
		enddo
		! swap line containing p with linej (incl. RHS)
!		call Mat_Swap(a,pline,j)
!		call Mat_Swap(inv,pline,j)
		
		! transform into upper triangular
		do i=j+1,nsize
			! (/((a(j,k)+1.3,j=1,1000),k=1,1000)/)
!			write(*,*) 'j=',j,'i=',i
			scalef = a(i,j) / a(j,j)
	!		write(*,*) 'scalef=',scalef
			do k=j,nsize
				if (i==2) then
!				write(*,*) 'a(k,i)',a(i,k),'- ',a(i,j),' / ',a(j,j),' * ',a(j,k)
				endif
				a(i,k) = a(i,k) - scalef * a(j,k) 
			enddo
			do k=1,nsize
				if (i==2) then
!				write(*,*) 'inv(i,k)::',inv(i,k),'- ',a(i,j),' / ',a(j,j),' * ',inv(j,k),'scalef=',scalef
				endif
				inv(i,k) = inv(i,k) - scalef * inv(j,k) 	
			enddo
		enddo
!		write(*,*) 'Matrix A after upper triangular, j=',j
!		call OutputMatrix(a)
!		write(*,*) 'Matrix INV after upper triangular, j=',j
!		call OutputMatrix(inv)
		! transform into lower triangular
		do i=j-1,1,-1
			! (/((a(j,k)+1.3,j=1,1000),k=1,1000)/)
			scalef = a(i,j) / a(j,j)
!			write(*,*) 'scalef=',scalef,' i=',i,'j=',j
			do k=j,nsize
!				write(*,*) 'a(k,i)',a(i,k),'- ',a(i,j),' / ',a(j,j),' * ',a(j,k)
				a(i,k) = a(i,k) - scalef * a(j,k) 
			enddo
			do k=1,nsize
!				write(*,*) 'inv(i,k)::',inv(i,k),'- ',a(i,j),' / ',a(j,j),' * ',inv(j,k),'scalef=',scalef
				inv(i,k) = inv(i,k) - scalef * inv(j,k) 	
			enddo
		enddo
!		write(*,*) 'Matrix A after lower triangular, j=',j
!		call OutputMatrix(a)
!		write(*,*) 'Matrix INV after lower triangular, j=',j
!		call OutputMatrix(inv)
		!normalize
		scalef = a(j,j)
		do k=1,nsize
	!		write(*,*) 'b(i,k)::',b(i,k),'- ',a(i,j),' / ',a(j,j),' * ',b(j,k),'scalef=',scalef
			inv(j,k) = inv(j,k) / scalef	
			a(j,k) = a(j,k) / scalef
!			write(*,*) 'a(j,k)',a(j,k),' / a(j,j)',a(j,j)
		enddo
	!	a(j,j) = a(j,j) / a(j,j)
	!	write(*,*) 'Matrix A NORMALIZED, j=',j
!		call OutputMatrix(inv)
	enddo
END SUBROUTINE Gauss_Inverse


!    Write a routine that performs a Gaussian Elimination with implicit pivoting to solve AX = B 
! for an arbitrary number of right-hand-side vectors b. The routine must read A from a file
! (as in the example given in section), and the rhs vectors B from another file. The solutions
!  must be printed to the screen by the calling program
! Returns the solutions in matrix B
SUBROUTINE Gauss_Elim(a,b,rowOrder)
	IMPLICIT NONE
	INTEGER :: n=0, pline=1, rhsVecCount,i,j,k
	real :: total,pmax,scalef
	real, dimension(:,:) :: a,b
	real, dimension(:), allocatable :: matScale
	real, dimension(:,:), allocatable :: rowOrder
	nsize = size(a,1)
	msize = size(a,2)
	rhsVecCount = size(b,2)
	if (size(b,1) .ne. nsize) then
		write(*,*) 'Arrays are not of equal size. Terminating.'
		stop
	endif
	allocate(matScale(nsize))
	allocate(rowOrder(nsize,1))
	total=1.
	do i=1,nsize
		rowOrder(i,1) = total
		total = total + 1.
	enddo
	
	write(*,*) '^ INPUT MATRICES ABOVE :: ','size(a,1)=',size(a,1),'size(b,1)=',size(b,1),'size(b,2)=',size(b,2)

	! (1) find scale factors for each row
	do i=1,nsize
		do j=1,msize
			if (matScale(i) < a(i,j)) then
				matScale(i) = a(i,j)
			endif
		enddo
	enddo
!	write(*,*) 'scale matrix:'
!	write (*,'(3f8.3)') matScale
	
	! (2) scale, pivot, and eliminate
	do j=1,nsize
		! find largest scaled pivot
		pmax=0
		pline=j
		do k=j,nsize
			if (pmax .lt. (a(k,j) / matScale(k))) then
				pmax = a(k,j) / matScale(k)
				pline = k
			endif
		enddo
		! swap line containing p with linej (incl. RHS)
		call Mat_Swap(rowOrder,pline,j)
		call Mat_Swap(a,pline,j)
		call Mat_Swap(b,pline,j)
		
		
		do i=j+1,nsize
			! (/((a(j,k)+1.3,j=1,1000),k=1,1000)/)
			scalef = a(i,j) / a(j,j)
	!		write(*,*) 'scalef=',scalef
			do k=j,nsize
	!			write(*,*) 'a(k,i)',a(i,k),'- ',a(i,j),' / ',a(j,j),' * ',a(j,k)
				a(i,k) = a(i,k) - scalef * a(j,k) 
			enddo
			do k=1,rhsVecCount
	!			write(*,*) 'b(i,k)::',b(i,k),'- ',a(i,j),' / ',a(j,j),' * ',b(j,k),'scalef=',scalef
				b(i,k) = b(i,k) - scalef * b(j,k) 	
			enddo
		enddo
		
	enddo
	
	! (3) backsubstitute
	
	do i=nsize,1,-1
		do k=1,rhsVecCount
			total=0
			do j=i+1,nsize
				total = total + (a(i,j) * b(j,k))
	!			write(*,*) 'total=',total
			enddo
			total = b(i,k) - total
			b(i,k) = total / a(i,i)
			
		enddo
	enddo

END SUBROUTINE Gauss_Elim


!    Write a routine that performs a Gaussian Elimination with implicit pivoting to solve AX = B 
! for an arbitrary number of right-hand-side vectors b. The routine must read A from a file
! (as in the example given in section), and the rhs vectors B from another file. The solutions
!  must be printed to the screen by the calling program
! Returns the solutions in matrix B
SUBROUTINE Mat_Swap(a,line1,line2)
	IMPLICIT NONE
	INTEGER :: line1, line2, nsize, i
	real :: temp
	real, dimension(:,:) :: a
	nsize = size(a,1)
	msize = size(a,2)
	do i=1,msize
		temp = a(line1,i)
		a(line1,i) = a(line2,i)
		a(line2,i) = temp
		!temp = a(i,line1)
		!a(i,line1) = a(i,line2)
		!a(i,line2) = temp
	enddo
END SUBROUTINE Mat_Swap


SUBROUTINE OutputMatrix(a)
	IMPLICIT NONE
	INTEGER :: i,j
	real, dimension(:,:) :: a
	do i=1,size(a,1)
		write(*,*) ( a(i,j), j=1,size(a,2) )
	enddo
END SUBROUTINE OutputMatrix

SUBROUTINE WriteMatrix(a,filename)
	IMPLICIT NONE
	INTEGER :: n,i,j
	CHARACTER(LEN=*) :: filename
	real, dimension(:,:) :: a
	n=18
	open(n,file=filename)
	write (n,*) size(a,1),size(a,2)
	do i=1,size(a,1)
		write(n,*) ( a(i,j), j=1,size(a,2) )
!		write(n,"100g15.5") ( a(i,j), j=1,size(a,2) )
	enddo
	!do i=1,size(a,1)
!		do j=1,size(a,2)
!			write(n,*) a(i,j)
	write(*,*) 'Printing output fit coordinates to: ',filename
	close(n)
END SUBROUTINE WriteMatrix

SUBROUTINE WriteVector(a,filename)
	IMPLICIT NONE
	INTEGER :: n,i,j
	CHARACTER(LEN=*) :: filename
	real, dimension(:) :: a
	n=18
	open(n,file=filename)
	write (n,*) size(a),'1'
	do i=1,size(a,1)
		write(n,*) a(i)
!		write(n,"100g15.5") ( a(i,j), j=1,size(a,2) )
	enddo
	!do i=1,size(a,1)
!		do j=1,size(a,2)
!			write(n,*) a(i,j)
	write(*,*) 'Printing output fit coordinates to: ',filename
	close(n)
END SUBROUTINE WriteVector

subroutine readandallocatemat(mat,filename)
	character (len=*) :: filename
	real, dimension(:,:), allocatable :: mat
	open(10,file=filename)
	read(10,*) nsize,msize
	allocate(mat(nsize,msize))
	do i=1,nsize
		read(10,*) ( mat(i,j), j=1,msize )
		write(*,*) ( mat(i,j), j=1,msize )
	enddo
	close(10)
	!write(*,*) 'Max value in matrix=',maxval(mat)
end subroutine readandallocatemat

real function trace(mat)
	real, dimension(nsize,msize) :: mat
	trace = 0.
	if(nsize.ne.msize) then
		write(*,*) 'This is not a square matrix, cannot calculate trace'
		else
		do i=1,nsize
			trace = trace + mat(i,i)
		enddo
	end if
end function trace


end module LinAl
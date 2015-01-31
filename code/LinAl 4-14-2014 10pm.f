module LinAl
implicit none
integer :: nsize,msize
integer,private :: i,j,k
contains

subroutine readandallocatemat(mat,filename)
	character*100 filename
	real, dimension(:,:), allocatable :: mat
	open(10,file=filename)
	read(10,*) nsize,msize
	allocate(mat(nsize,msize))
	do i=1,nsize
		read(10,*) ( mat(i,j), j=1,msize )
		write(*,*) ( mat(i,j), j=1,msize )
	enddo
	close(10)
	write(*,*) 'Max value in matrix=',maxval(mat)
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

  subroutine inverse(a,inv)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
		IMPLICIT NONE
	INTEGER :: n,i,j,k
	real :: total,pmax,scalef
	real, dimension(:,:) :: a,L,U
	real, dimension(:,:), allocatable :: inv
	nsize = size(a,1)
	msize = size(a,2)

	allocate(inv(nsize,msize))
	allocate(L(nsize,msize))
	allocate(U(nsize,msize))
	
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
	

	! step 0: initialization for matrices L and U and b
	! Fortran 90/95 aloows such operations on matrices
	L=0.0
	U=0.0
	b=0.0

	! step 1: forward elimination
	do k=1, n-1
	   do i=k+1,n
		  coeff=a(i,k)/a(k,k)
		  L(i,k) = coeff
		  do j=k+1,n
			 a(i,j) = a(i,j)-coeff*a(k,j)
		  end do
	   end do
	end do

	! Step 2: prepare L and U matrices 
	! L matrix is a matrix of the elimination coefficient
	! + the diagonal elements are 1.0
	do i=1,n
	  L(i,i) = 1.0
	end do
	! U matrix is the upper triangular part of A
	do j=1,n
	  do i=1,j
		U(i,j) = a(i,j)
	  end do
	end do

	! Step 3: compute columns of the inverse matrix C
	do k=1,n
	  b(k)=1.0
	  d(1) = b(1)
	! Step 3a: Solve Ld=b using the forward substitution
	  do i=2,n
		d(i)=b(i)
		do j=1,i-1
		  d(i) = d(i) - L(i,j)*d(j)
		end do
	  end do
	! Step 3b: Solve Ux=d using the back substitution
	  x(n)=d(n)/U(n,n)
	  do i = n-1,1,-1
		x(i) = d(i)
		do j=n,i+1,-1
		  x(i)=x(i)-U(i,j)*x(j)
		end do
		x(i) = x(i)/u(i,i)
	  end do
	! Step 3c: fill the solutions x(n) into column k of C
	  do i=1,n
		c(i,k) = x(i)
	  end do
	  b(k)=0.0
	end do
end subroutine inverse

!    Write a routine that performs an LU decomposition, and returns the LU form of the matrix. 
!     The LU form must be printed to the screen by the calling program. Pivoting is not necessary.
!    (Optional, could be useful for you later): Copy the LU decomposition and backsubstitution routines from
!     Numerical Recipes or any other source, and learn how to use it to solve linear algebra problems. 
!     You're free to choose a F77 or an F90 source.
SUBROUTINE LU_Decomposition(a,L,U)
	IMPLICIT NONE
	INTEGER :: n=0, pline=1, rhsVecCount,i,j,k
	real :: total,pmax,scalef
	real, dimension(:,:) :: a
	real, dimension(:,:), allocatable :: L,U
	nsize = size(a,1)
	msize = size(a,2)
	allocate(L(nsize,msize))
	allocate(U(nsize,msize))
	
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
	write(*,*) 'scale matrix:'
	write (*,'(3f8.3)') matScale
	
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
	!	call Mat_Swap(rowOrder,pline,j)
		call Mat_Swap(a,pline,j)
		call Mat_Swap(inv,pline,j)
		
		! transform into upper triangular
		do i=j+1,nsize
			! (/((a(j,k)+1.3,j=1,1000),k=1,1000)/)
			scalef = a(i,j) / a(j,j)
	!		write(*,*) 'scalef=',scalef
			do k=j,nsize
	!			write(*,*) 'a(k,i)',a(i,k),'- ',a(i,j),' / ',a(j,j),' * ',a(j,k)
				a(i,k) = a(i,k) - scalef * a(j,k) 
			enddo
			do k=1,nsize
	!			write(*,*) 'b(i,k)::',b(i,k),'- ',a(i,j),' / ',a(j,j),' * ',b(j,k),'scalef=',scalef
				inv(i,k) = inv(i,k) - scalef * inv(j,k) 	
			enddo
		enddo
		
		! transform into lower triangular
		do i=j-1,1,-1
			! (/((a(j,k)+1.3,j=1,1000),k=1,1000)/)
			scalef = a(i,j) / a(j,j)
			write(*,*) 'scalef=',scalef,' i=',i,'j=',j
			do k=j,nsize
				write(*,*) 'a(k,i)',a(i,k),'- ',a(i,j),' / ',a(j,j),' * ',a(j,k)
				a(i,k) = a(i,k) - scalef * a(j,k) 
			enddo
			do k=1,nsize
	!			write(*,*) 'b(i,k)::',b(i,k),'- ',a(i,j),' / ',a(j,j),' * ',b(j,k),'scalef=',scalef
				inv(i,k) = inv(i,k) - scalef * inv(j,k) 	
			enddo
		enddo
		
		!normalize
		do i=1,nsize
	!		write(*,*) 'b(i,k)::',b(i,k),'- ',a(i,j),' / ',a(j,j),' * ',b(j,k),'scalef=',scalef
			inv(i,j) = inv(i,j) / a(j,j)	
		enddo
		a(j,j) = a(j,j) / a(j,j)
		
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
	
	write(*,*) '^ INPUT MATRICES ABOVE :: ','size(a,1)=',size(a,1),'size(b,1)=',size(b,1)

	! (1) find scale factors for each row
	do i=1,nsize
		do j=1,msize
			if (matScale(i) < a(i,j)) then
				matScale(i) = a(i,j)
			endif
		enddo
	enddo
	write(*,*) 'scale matrix:'
	write (*,'(3f8.3)') matScale
	
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

end module LinAl
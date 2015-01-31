!  *********! SECTION 3 --------- AMS 213 **************! !

program section3
use LinAl
use NumRec
implicit none
real,dimension(:,:), allocatable :: matA,matB,matRowOrder,matInv,matL,matU,matLU!,matScale
real,dimension(:), allocatable :: matPivot

character*100 filename
integer :: i,j,k,n,m,info
real :: x=1.,y,guess=1.

if(iargc().ne.2) then
	write(*,*) 'Wrong number of arguments (need file names for LHS (A) and RHS (B)'
	stop
endif


! Gaussian Elimination w/ Partial Pivoting
write(*,*) '*******! GaussJ.f90 : Numerical Recipes *******! '
! Allocate matrices A,B, and scale
call getarg(1,filename)
write(*,*) '*******! Input Matrix A *******! '
call readandallocatemat(matA,filename)
call getarg(2,filename)
write(*,*) '*******! Input Matrix B *******! '
call readandallocatemat(matB,filename)

call gaussj(matA,matB)
write(*,*) '*******! Output Matrix A *******! '
call OutputMatrix(matA)
write(*,*) '*******!  Output Matrix B *******! '
call OutputMatrix(matB)

! Gaussian Elimination w/ Partial Pivoting
write(*,*) '*******! GaussJ.f90 : Numerical Recipes *******! '
deallocate(matA)
deallocate(matB)
! Allocate matrices A,B, and scale
call getarg(1,filename)
write(*,*) '*******! Input Matrix A *******! '
call readandallocatemat(matA,filename)
call getarg(2,filename)
write(*,*) '*******! Input Matrix B *******! '
call readandallocatemat(matB,filename)

call gaussj(matA(1:size(matA,1),1:size(matA,2)), matB(1:size(matB,1),1:size(matB,2)))
write(*,*) '*******! Output Matrix A *******! '
call OutputMatrix(matA)
write(*,*) '*******!  Output Matrix B *******! '
call OutputMatrix(matB)

! Matrix inverse calculation
write(*,*) '*******! Matrix inverse calculation *******! '
! Allocate matrices A
call getarg(1,filename)
deallocate(matA)
write(*,*) '*******! Input Matrix A *******! '
call readandallocatemat(matA,filename)

call Gauss_Inverse(matA,matInv)
write(*,*) '*******!  Output Matrix A *******! ' 
call OutputMatrix(matA)
write(*,*) '*******!  Output Matrix A Inverse *******! '
call OutputMatrix(matInv)
call WriteMatrix(matInv,'matInv.dat')


! Matrix inverse calculation
write(*,*) '*******! Matrix inverse calculation *******! '
! Allocate matrices A
call getarg(1,filename)
deallocate(matA)
write(*,*) '*******! Input Matrix A *******! '
call readandallocatemat(matA,filename)

!  SGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P ! L ! U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.

allocate matPivot(size(matA,1))
SUBROUTINE SGETRF( size(matA,1), size(matA,2), matA, size(matA,1), matPivot, info )


!  SGETRI computes the inverse of a matrix using the LU factorization
!  computed by SGETRF.
!  This method inverts U and then computes inv(A) by solving the system
!  inv(A)*L = inv(U) for inv(A).
!  Arguments
!  =========
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the factors L and U from the factorization
!          A = P*L*U as computed by SGETRF.
!          On exit, if INFO = 0, the inverse of the original matrix A.
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from SGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
!          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimal performance LWORK >= N*NB, where NB is
!          the optimal blocksize returned by ILAENV.
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!                singular and its inverse could not be computed.


!SUBROUTINE SGETRI( 3, matA, 3, IPIV, WORK, LWORK, info )
write(*,*) '*******!  Output Matrix A *******! ' 
call OutputMatrix(matA)
write(*,*) '*******!  Output Matrix A Inverse *******! '
call OutputMatrix(matInv)
call WriteMatrix(matInv,'matInv.dat')






end program section3


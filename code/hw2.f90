!  ********** HOMEWORK 2 --------- AMS 213 *************** !

program hw2
use LinAl
implicit none
real,dimension(:,:), allocatable :: matA,matB,matRowOrder,matInv,matL,matU,matLU!,matScale
character*100 filename
!integer, parameter :: dimmat = 3, rhsVecCount = 1
!real, dimension(dimmat,dimmat) :: a
!real, dimension(dimmat,rhsVecCount) :: b
!real, dimension(dimmat) :: matScale
!real, dimension(:), allocatable :: matScale
integer :: i,j,k
real :: x=1.,y,guess=1.

if(iargc().ne.2) then
	write(*,*) 'Wrong number of arguments (need file names for LHS (A) and RHS (B)'
	stop
endif


! Gaussian Elimination w/ Partial Pivoting
write(*,*) '******** Gaussian Elimination w/ Partial Pivoting ******** '
! Allocate matrices A,B, and scale
call getarg(1,filename)
write(*,*) '******** Input Matrix A ******** '
call readandallocatemat(matA,filename)
call getarg(2,filename)
write(*,*) '******** Input Matrix B ******** '
call readandallocatemat(matB,filename)

call Gauss_Elim(matA,matB,matRowOrder)
write(*,*) '******** Output Matrix A ******** '
call OutputMatrix(matA)
write(*,*) '********  Output Matrix B ******** '
call OutputMatrix(matB)
write(*,*) '********  Output Row Ordering ******** '
call OutputMatrix(matRowOrder)

! Matrix inverse calculation
write(*,*) '******** Matrix inverse calculation ******** '
! Allocate matrices A
call getarg(1,filename)
deallocate(matA)
write(*,*) '******** Input Matrix A ******** '
call readandallocatemat(matA,filename)

call Gauss_Inverse(matA,matInv)
write(*,*) '********  Output Matrix A ******** ' 
call OutputMatrix(matA)
write(*,*) '********  Output Matrix A Inverse ******** '
call OutputMatrix(matInv)
call WriteMatrix(matInv,'matInv.dat')

! LU Decomposition
write(*,*) '******** LU Decomposition ******** '
write(*,*) '******** Input Matrix A ******** '
deallocate(matA)
! Allocate matrices A,B, and scale
call getarg(1,filename)
call readandallocatemat(matA,filename)
call LU_Decomposition(matA,matL,matU)
write(*,*) '********  Output Matrix L ******** '
call OutputMatrix(matL)
write(*,*) '********  Output Matrix U ******** '
call OutputMatrix(matU)
matLU = matmul(matL,matU)
write (*,*) '********  Matrix L*U :: matmul(matL,matU) ******** '
call OutputMatrix(matLU)
!call Mat_Swap(matA,1,2)
!write (*,'(3f8.3)') matA
!write (*,'(3f8.3)') matB

end program hw2


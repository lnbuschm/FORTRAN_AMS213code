!  ********** HOMEWORK 3 --------- AMS 213 *************** !

program hw3
use LinAl
implicit none
real,dimension(:,:), allocatable :: matA,matB,matC,matRowOrder,matInv,matL,matU,matLU
real :: ax,bx,cx
character*100 filename,filename2
integer :: i,j,k
real :: x=1.,y,guess=1.

write(*,*) '******** HOMEWORK 3  --  AMS213 ******** '

! Cholesky Decomposition
write(*,*) '******** Cholesky Decomposition ******** '

! Allocate matrices A,B, and scale
call getarg(1,filename)
write(*,*) '******** Input Matrix A ******** '
call readandallocatemat(matA,filename)
call getarg(2,filename)
write(*,*) '******** Input Matrix B ******** '
call readandallocatemat(matB,filename)

! Solve Cholesky decomposition of A, store in C
call CHOLESKYD(matA,matC)
write(*,*) '********  Output Matrix C : Cholesky Decomposition of A ******** '
call OutputMatrix(matC)
deallocate(matC)

! Solve Ax=B
call CHOLESKYDSOL(matA,matB)
write(*,*) '********  Output Matrix B : Solution to Ax=B ******** '
call OutputMatrix(matB)
 
call getarg(3,filename)
write(*,*) '******** Cubic Least Squares ******** '

!call CHOLESKYLEASTSQUARES(ax,bx,cx,3,filename,'atkinsonfit.dat')

call CUBICLEASTSQUARES(ax,bx,cx,filename,'atkinsonfit.dat')

call LAPACKLEASTSQUARES(ax,bx,cx,filename,'atkinsonfitLAPACK.dat')


end program hw3


!  ********** MIDTERM PROJECT --------- AMS 213 *************** !

program project1
use LinAl
implicit none
real,dimension(:,:), allocatable :: matA,matB,matC,matRowOrder,matInv,matL,matU,matLU,matQ,matR
real,dimension(:), allocatable :: c,d,parameters
real :: ax,bx,cx
character (len=100) :: filename
integer :: i,j,k,order=6
real :: x=1.,y,guess=1.
LOGICAL :: sing

write(*,*) '******** MIDTERM PROJECT  --  AMS213 ******** '

! Allocate matrices A,B, and scale
call getarg(1,filename)
write(*,*) '******** Input Matrix A ******** '
call readandallocatemat(matA,filename)


! Test canned QR Decomposition routine 
!allocate(c(nsize))
!allocate(d(nsize))
!call qrdcmp(matA,size(matA,1),size(matA,1),c,d,sing)
!write(*,*) 'Vector C',c
!write(*,*) 'Vector D',d
!call OutputMatrix(matA)
!deallocate(matA)
!call getarg(1,filename)
!write(*,*) '******** Input Matrix A ******** '
!call readandallocatemat(matA,filename)
!deallocate(c)
!deallocate(d)

! Solve QR Decomposition of A
call QRDECOMP(matA,matQ,matR)  
write(*,*) '********  Output Matrix Q : QR Decomposition of A ******** '
call OutputMatrix(matQ)
write(*,*) '********  Output Matrix R : QR Decomposition of A ******** '
call OutputMatrix(matR)

!call QRDECOMP(matA,c,d)
!write(*,*) '********  Output Matrix A : QR Decomposition of A ******** '
!call OutputMatrix(matA)

! Solve Ax=B using QR Decomposition & Backsubstitution
call getarg(2,filename)
write(*,*) '******** Input Matrix B ******** '
call readandallocatemat(matB,filename)
call QRSOLVE(matA,matB)
write(*,*) '******** Output Solution Matrix ******** '
call OutputMatrix(matB)


write(*,*) '********  QR LEAST SQUARES ******** '
call QRLEASTSQUARES(order,parameters,'atkinson.dat','qrparameters.dat','qrfit.dat')

deallocate(parameters)
write(*,*) '********  CHOLESKY LEAST SQUARES ******** '
call CHOLESKYLEASTSQUARES(order,parameters,'atkinson.dat','choleskyparameters.dat','choleskyfit.dat')


! Solve Cholesky decomposition of A, store in C
!call CHOLESKYD(matA,matC)
!write(*,*) '********  Output Matrix C : Cholesky Decomposition of A ******** '
!call OutputMatrix(matC)
!deallocate(matC)

! Solve Ax=B
!call CHOLESKYDSOL(matA,matB)
!write(*,*) '********  Output Matrix B : Solution to Ax=B ******** '
!call OutputMatrix(matB)
 
call getarg(3,filename)
!write(*,*) '******** Cubic Least Squares ******** '
!call CUBICLEASTSQUARES(ax,bx,cx,filename,'atkinsonfit.dat')

!call LAPACKLEASTSQUARES(ax,bx,cx,'atkinson.dat','atkinsonfitLAPACK.dat')


end program project1


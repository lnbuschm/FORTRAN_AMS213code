!  ********** HOMEWORK 1 --------- AMS 213 *************** !
!   Write a program that calls this function and prints the result to the screen

program hw1
implicit none
real :: x=1.,y,guess=1.
real :: fcosx,fcosprimex,numPrecision,func1,func1prime
external fcosx,fcosprimex,func1,func1prime

!write(*,*) 'Numerical Precision',numPrecision() 
x=numPrecision()



call Newton_Raphson(y,fcosx,fcosprimex,guess)
write(*,*) 'Subroutine:: cos(x)=0, x=',y,'= pi/2'

call Newton_Raphson(y,func1,func1prime,guess)
write(*,*) 'Subroutine:: exp(x)-2=0, x=',y


!call Newton_Raphson(cos,0.0)
!call NR(sol,fcosx,guess)  !  where sol is solution, fcosx is function written in tutorial, guess is =1
!write(*,*)
!  subroutine fcosx(x,f,fprime  !takes in x, returns f and fprime ,,, f=cosx, fprime=-sinx
!or call NR(Sol,fcosx,fprime,guess)

end program hw1



! Write a function that returns the precision of the machine you are using for single and double-precision real. 
! Hint: use the fact that 1+x only returns something different from 1 if x is larger than the numerical precision.
function numPrecision()
implicit none
double precision :: y=1.
real :: x=1,numPrecision
!real numPrecision
do while (x+1 /= 1)
    x=x/2.
!	if (x+1. /= 1.) write(*,*) x+1.,' /= 1'
enddo
write(*,*)  'Numerical Precision: real',x
do while (y+1 /= 1)
    y=y/2.
!	if (x+1. /= 1.) write(*,*) x+1.,' /= 1'
enddo
write(*,*)  'Numerical Precision: double',y
!return x
numPrecision = x;
end function numPrecision

function fcosx(x)
implicit none
real :: fcosx,x
fcosx = cos(x)
end function fcosx

function fcosprimex(x)
implicit none
real :: fcosprimex,x
fcosprimex = -1 * sin(x)
end function fcosprimex

function func1(x)
implicit none
real :: func1,x
func1 = exp(x)-2.
end function func1

function func1prime(x)
implicit none
real :: func1prime,x
func1prime = exp(x)
end function func1prime

!-----Newton-Raphson---------------------------------------------------
!   Write a program that performs a Newton-Raphson search for the solution of the equation cos(x) = 0 
! to machine precision starting from the guess x=1. You are free to choose to use either single (real)
! or double (real*8) precision. The Newton-Raphson algorithm must be written in a separate routine, 
! and so is the function cos(x). The calling program should call the Newton-Raphson solver and pass 
! the function cos(x) as an argument then write the solution to the screen. You MUST write your own routine,
! although you may inspire yourself from the one given in Numerical Recipes or other computing course.
!---------------------------------------------------------------------
!call NR(sol,fcosx,guess)  !  where sol is solution, fcosx is function written in tutorial, guess is =1

SUBROUTINE Newton_Raphson(sol,func,funcprime,guess)
IMPLICIT NONE
INTEGER :: n
REAL :: guess,func,funcprime,xn=0
REAL, INTENT(OUT) :: sol
external func,funcprime
!real func
sol = guess
n=0
do while (n .lt. 100)
!write(*,*) 'sol=',sol
    sol = sol - func(sol)/funcprime(sol)
    n=n + 1
!	if (x+1. /= 1.) write(*,*) x+1.,' /= 1'
enddo
!write(*,*) 'Subroutine:: cos(x)=0, x=',sol,'= pi/2'
!Area = Pi * r * r

END SUBROUTINE Newton_Raphson
program nameofprogram
implicit none
integer :: i,j,k,n
real :: x,y,z
n=10
x = 3.61
y = cos(x)
z = x + y
i = 3
j = i**2
k = i - j
open(n,file='mydata.dat')
write(n,*) i
write(n,*) j
write(n,*) k
write(n,*) x
write(n,*) y
write(n,*) z
close(n)
end program nameofprogram
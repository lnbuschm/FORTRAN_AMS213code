program nameofprogram
implicit none
integer :: i,j,k
real :: x,y,z
x = 3.61
y = cos(x)
z = x + y
i = 3
j = i**2
k = i - j
write(*,*) 'The value of x is ',x,' and the value of y is', y
end program nameofprogram
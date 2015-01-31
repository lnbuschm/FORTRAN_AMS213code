program firstLinAlprog
use LinAl
real,dimension(:,:), allocatable :: mat
real :: x
character*100 filename
if(iargc().ne.1) then
write(*,*) 'Wrong number of arguments (need file name)'
stop
endif
call getarg(1,filename)
call readandallocatemat(mat,filename)
x = trace(mat)
write(*,*) 'The trace of this matrix is ', x
deallocate(mat)
end program firstLinAlprog

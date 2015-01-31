silly:
	gfortran section2a.f90 silly
hw1:
	gfortran hw1.f90 -o hw1
hw2:
	gfortran LinAl.f90 hw2.f90 -o hw2
	./hw2 Amat.dat Bmat.dat
hw3:
	gfortran LinAl.f90  hw3.f90  -lliblapack -llibblas -o hw3  
	./hw3 hw3_a.dat hw2_b.dat atkinson.dat
hw4:
	gfortran LinAl.f90 DiffEq.f90 hw4.f90  -lliblapack -llibblas -o hw4
	./hw4
hw5:
	gfortran LinAl.f90 DiffEq.f90 hw5.f90  -lliblapack -llibblas -o hw5
	./hw5
hw6:
	gfortran LinAl.f90 DiffEq.f90 hw6.f90  -lliblapack -llibblas -o hw6
	./hw6
projectfinal:
	gfortran LinAl.f90 DiffEq.f90 projectfinal.f90  -lliblapack -llibblas -o projectfinal
	./projectfinal
amsteam:
	gfortran teammod.f90 deriv.f90 teamdriver.f90  -o amsteam
	./amsteam
project1:
	gfortran LinAl.f90 project1.f90  -lliblapack -llibblas -o project1 
	./project1 hw3_a.dat hw2_b.dat atkinson.dat
section3:
	gfortran LinAl.f90 gaussj.f90 NumRec.f90 nrtype.f90 nrutil.f90 section3.f90 -o section3
	./section3 hw2_a.dat hw2_b.dat
section3a:
	gfortran LinAl.f90 gaussj.f90 NumRec.f90 nrtype.f90 nrutil.f90 section3.f90 -lliblapack -llibblas  -o section3a
	./section3a hw2_a.dat hw2_b.dat
lapacktest:
	gfortran lapacktest.f90 -lliblapack -llibblas -o lapacktest
.PHONY: hw1 hw2 hw3 lapacktest section3 section3a project1 hw4 hw5 hw6 projectfinal
module NumRec

contains

SUBROUTINE gaussj(a,b)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
	INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
	LOGICAL(LGT), DIMENSION(size(a,1)) :: lpiv
	REAL(SP) :: pivinv
	REAL(SP), DIMENSION(size(a,1)) :: dumc
	INTEGER(I4B), TARGET :: irc(2)
	INTEGER(I4B) :: i,l,n
	INTEGER(I4B), POINTER :: irow,icol
	n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
	irow => irc(1)
	icol => irc(2)
	ipiv=0
	do i=1,n
		lpiv = (ipiv == 0)
		irc=maxloc(abs(a),outerand(lpiv,lpiv))
		ipiv(icol)=ipiv(icol)+1
		if (ipiv(icol) > 1) call nrerror('gaussj: singular matrix (1)')
		if (irow /= icol) then
			call swap(a(irow,:),a(icol,:))
			call swap(b(irow,:),b(icol,:))
		end if
		indxr(i)=irow
		indxc(i)=icol
		if (a(icol,icol) == 0.0) &
			call nrerror('gaussj: singular matrix (2)')
		pivinv=1.0_sp/a(icol,icol)
		a(icol,icol)=1.0
		a(icol,:)=a(icol,:)*pivinv
		b(icol,:)=b(icol,:)*pivinv
		dumc=a(:,icol)
		a(:,icol)=0.0
		a(icol,icol)=pivinv
		a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
		b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
		a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
		b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
	end do
	do l=n,1,-1
		call swap(a(:,indxr(l)),a(:,indxc(l)))
	end do
	END SUBROUTINE gaussj

SUBROUTINE ludcmp(a,indx,d)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(SP), INTENT(OUT) :: d
	REAL(SP), DIMENSION(size(a,1)) :: vv
	REAL(SP), PARAMETER :: TINY=1.0e-20_sp
	INTEGER(I4B) :: j,n,imax
	n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
	d=1.0
	vv=maxval(abs(a),dim=2)
	if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
	vv=1.0_sp/vv
	do j=1,n
		imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
		if (j /= imax) then
			call swap(a(imax,:),a(j,:))
			d=-d
			vv(imax)=vv(j)
		end if
		indx(j)=imax
		if (a(j,j) == 0.0) a(j,j)=TINY
		a(j+1:n,j)=a(j+1:n,j)/a(j,j)
		a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
	end do
	END SUBROUTINE ludcmp

	SUBROUTINE lubksb(a,indx,b)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,n,ii,ll
	REAL(SP) :: summ
	n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
	ii=0
	do i=1,n
		ll=indx(i)
		summ=b(ll)
		b(ll)=b(i)
		if (ii /= 0) then
			summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
		else if (summ /= 0.0) then
			ii=i
		end if
		b(i)=summ
	end do
	do i=n,1,-1
		b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
	end do
	END SUBROUTINE lubksb

SUBROUTINE choldc(a,p)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: p
	INTEGER(I4B) :: i,n
	REAL(SP) :: summ
	n=assert_eq(size(a,1),size(a,2),size(p),'choldc')
	do i=1,n
		summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
		if (summ <= 0.0) call nrerror('choldc failed')
		p(i)=sqrt(summ)
		a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
	end do
	END SUBROUTINE choldc

	SUBROUTINE cholsl(a,p,b,x)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(IN) :: p,b
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	INTEGER(I4B) :: i,n
	n=assert_eq((/size(a,1),size(a,2),size(p),size(b),size(x)/),'cholsl')
	do i=1,n
		x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/p(i)
	end do
	do i=n,1,-1
		x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/p(i)
	end do
	END SUBROUTINE cholsl


end module NumRec

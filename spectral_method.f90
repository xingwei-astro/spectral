!!! slove equation d^u/dy^2=0 with fixed boundary conditions, x is in [-1,1], y=(x+1)/2 is in [0,1]
!!! ansatz u=Simga_0^(n+1) (uhat_j*T_j(x_i)),  i, j=0, ..., n+1
!!! gfortran -o laplace spectral_method.f90 solve.f90

module globe
implicit none
double precision pi
parameter (pi=acos(-1.d0))
integer n
parameter (n=10)
double precision x(0:n+1), y(0:n+1), u(0:n+1)   ! physical space 
double precision uhat(0:n+1)   ! spectral coefficients
double precision u1, u2   ! boundary conditions
parameter (u1=1.d0, u2=2.d0)
double precision a(0:n+1,0:n+1)   ! left-hand-side matrix
double precision b(0:n+1)   ! right-hand-side vector
end module globe

program main
use globe
implicit none
integer i, j   ! i is index of collocation points, j is index of spectral coefficients
double precision TTT
double precision a_inv(0:n+1,0:n+1)

x(0)=-1.d0
x(n+1)=1.d0
do i=1, n
 x(i)=-cos(dfloat(2*i-1)*pi/dfloat(2*n))   ! inner collocation points, zeros of T_n
! write(6,*) x(i), TTT(0,n,x(i))
enddo
do i=0, n+1
 y(i)=(x(i)+1.d0)/2.d0
! write(6,*) y(i)
enddo

! collocate equation on inner points
do j=0, n+1
 do i=1, n
  a(i,j)=TTT(2,j,x(i))
 enddo
enddo
! collocate equation on boundary points
do j=0, n+1
 a(0,j)=TTT(0,j,x(0))
 a(n+1,j)=TTT(0,j,x(n+1))
enddo
! the right-hand-side vector
do i=1, n
 b(i)=0.d0
enddo
b(0)=u1
b(n+1)=u2

! solve a*uhat=b
call r8mat_fs(n+2,a,b,uhat)

! calculate u
do i=0, n+1
 u(i)=0.d0
 do j=0, n+1
  u(i)=u(i)+uhat(j)*TTT(0,j,x(i))
 enddo
enddo

! output
open(1,file='laplace.dat',form='formatted')
do i=0, n+1
 write(1,'(2F10.3)') y(i), u(i)
enddo
close(1)
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!   TTT(K,M,X) = the K-th derivative of Tm(X), the M-th Chebyshev polynomial evaluated at X.

      DOUBLE PRECISION FUNCTION TTT(K,M,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(0:10000), B(0:10000)
      DO 10 I=0,M-1
       A(I)=0.D0
10    CONTINUE
      A(M)=1.D0
      DO 20 J=1,K
       CALL DIFF(A,M)
20    CONTINUE
      B(M+2)=0.D0
      B(M+1)=0.D0
      DO 30 I=M,0,-1
       B(I)=2.D0 * X * B(I+1)  -  B(I+2)  +  A(I)
30    CONTINUE
      TTT=(A(0) + B(0) - B(2))/2.D0
      RETURN
      END
      SUBROUTINE DIFF(A,M)
      DOUBLE PRECISION A(0:10000), C(0:10000)
      C(M+1)=0.D0
      C(M)=0.D0
      DO 10 I=M-1,0,-1
       C(I)=C(I+2)  +  DFLOAT(2*I+2) * A(I+1)
10    CONTINUE
      C(0)=C(0)/2.D0
      DO 20 I=0,M
       A(I)=C(I)
20    CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
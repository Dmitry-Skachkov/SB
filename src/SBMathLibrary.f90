


 Module SBMathLibrary                                               ! standard subroutines for integration and splines
      implicit none
 Contains




      SUBROUTINE QSL3D(RINT,A,B,F,E)
!
!     Adaptive integration using Lobatto formula of 11th power with guaranteed accuracy. 
!
!     RINT - on exit contains the integral value
!     A, B - lower and upper limits for integration 
!     F    - external subroutine with integration function
!     E    - accuracy of calculation
!     In a case when the integration function has singularities and the necessary accuracy can not be achieved 
!     the program may stop with the corresponding message.
!
! Copyright (C) 1974-2021 Library of the numerical analysis, Research Computing Center, 
!                         Lomonosov Moscow State University
!                         http://num-anal.srcc.msu.ru/lib_na/libnal.htm
!     The library can be used without any limitations only for non-commercial usage.
!     If the Library or any part of it is installed in organization for use in a shared access mode, 
!     then the organization must arrange non-commercial usage of the Library's subprograms.
!     Users or organizations have no right to transfer the programs of the Library or their modifications 
!     to other persons or organizations for non-commercial usage without written permission 
!     from the Research Computing Center of Moscow State University.
!     For commercial usage please contact Research Computing Center at <arush@srcc.msu.ru>
!
      External F
      DIMENSION AY(7),BY(7),KSI(7),OM(3),D(3)
      INTEGER OM,N,IL,O
      real(8) RINT,A,B,F,E,XA,AY,BY,KSI,D,XAA,INT,INT2,I,J,K,L,M,C,W,E1,Z2,H,Z1,EPS,P,H0,Q,R,SYS051
      LOGICAL NI
      DATA AY/0.023809523809D0,0.243809523809D0,0.023809523809D0,0.138413023681D0,0.215872690605D0,0.215872690605D0,0.138413023681D0/
      DATA BY/-0.076190474D0,0.243809523809D0,-0.076190474D0,0.183701962D0,-0.229416249D0,-0.229416249D0,0.183701962D0/
      DATA KSI/0.D0,0.5D0,1.D0,0.084888051861D0,0.265575603265D0,0.734424396735D0,0.915111948139D0/
      DATA SYS051/1.1107652D-16/
      EPS=E
      IL=1
      NI=A.GT.B
      IF(.NOT.(NI)) GO TO 1
      C=A
      A=B
      B=C
    1 W=B-A
      IF(.NOT.(W.EQ.0.D0)) GO TO 2
      INT=0.D0
      E1=0.D0
      GO TO 28
    2 Z2=E/W
      H=1.D0
      N=0
      Z1=EPS/W
    3 P=0.D0
      E1=0.D0
      INT=0.D0
      OM(1)=0
      OM(2)=0
      OM(3)=0
    4 K=P+H
      L=0.D0
      M=0.D0
      DO 7 O=1,3
      IF(OM(O).EQ.0) GO TO 5
      I=D(O)
      GO TO 6
    5 J=H*KSI(O)+P
      D(O)=F(W*J+A)
      I=D(O)
      OM(O)=1
    6 L=L+AY(O)*I
      M=M+BY(O)*I
    7 CONTINUE
      DO 8 O=4,7
      J=H*KSI(O)+P
      I=F(W*J+A)
      L=L+AY(O)*I
      M=M+BY(O)*I
    8 CONTINUE
      L=H*L
      M=DABS(H*M)
      IF(M-Z1) 9,18,18
    9 P=K
      E1=E1+M
      INT=INT+L
      I=1.D0-P
      N=N+1
      IF(N-1) 11,10,11
   10 H0=H
   11 IF(1.D0+I/64.D0-1.D0) 12,12,14
   12 J=DABS(F(B))+DABS(D(3))
      J=J*I
      IF(J-Z2) 22,22,13
   13 XA=A+P*W
      INT=INT*W
      XAA=XA
      print 111,XAA
111   format('The accuracy can not be achived at singular point',E20.12)
      GO TO 28
   14 D(1)=D(3)
      OM(2)=0
      OM(3)=0
      OM(1)=1
      IF(M-Z1/2.D0**7) 15,16,16
   15 H=2.D0*H
   16 IF(H-I) 4,4,17
   17 H=I
      GO TO 4
   18 H=H/2.D0
      IF(H-SYS051*P) 19,19,21
   19 IF(P) 20,20,13
   20 IF(H/4.D0) 13,13,21
   21 D(3)=D(2)
      OM(2)=0
      GO TO 4
   22 Q=Z2*(1.D0+DABS(INT))
      IF(E1-Q) 27,27,23
   23 IF(IL-1) 24,25,24
   24 R=DABS(INT-INT2)
      IF(R-Q) 26,26,25
   25 INT2=INT
      IL=2
      EPS=EPS/2.D0**7
      H=H0
      GO TO 3
   26 E1=R
   27 INT=INT*W
      XA=B
      IF(.NOT.(NI)) GO TO 28
      C=B
      B=A
      A=C
      INT=-INT
   28 RINT=INT
      RETURN
      END SUBROUTINE QSL3D





   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  
!  Author: Alex Godunov, January 2010 
!  https://ww2.odu.edu/~agodunov/computing/programs/
!  spline.f90 program is based on http://netlib.org/ fortran version of program spline.f
!  the accompanying function ispline can be used for interpolation
!======================================================================
integer     :: n
real(8)     :: x(n), y(n), b(n), c(n), d(n)
integer     :: i, j, gap
real(8)     :: h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination 
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline




  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!
!  Author: Alex Godunov, January 2010 
!  https://ww2.odu.edu/~agodunov/computing/programs/
!=======================================================================
real(8)        :: ispline
integer        :: n
real(8)        :: u, x(n), y(n), b(n), c(n), d(n)
integer        :: i, j, k
real(8)        :: dx

!*** modified by Dmitry Skachkov: adding extrapolation for points outside of X(1)..X(N) 
! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
!  ispline = y(1)
  dx = u - x(1)                                          !!***
  ispline = y(1) + dx*(b(1) + dx*(c(1) + dx*d(1)))       !!***
  return
end if
if(u >= x(n)) then
!  ispline = y(n)
  dx = u - x(n)                                          !!***
  ispline = y(n) + dx*(b(n) + dx*(c(n) + dx*d(n)))       !!***
  return
end if
!***

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline





 end Module SBMathLibrary







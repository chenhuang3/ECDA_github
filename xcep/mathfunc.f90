module mathfunc

contains

function h00(n1,n2,n3,t)
!---------------------------------------------------
! hermite function basis h00
! Used for cubic Hermite spline
!---------------------------------------------------
  IMPLICIT NONE

  integer, parameter :: dp=8
  integer :: n1,n2,n3
  REAL(kind=dp),intent(in)  :: t(n1,n2,n3)
  REAL(kind=dp)             :: h00(n1,n2,n3)

  h00 = 2._DP * t**3 - 3._DP*t**2 + 1._DP
  return
END FUNCTION h00

FUNCTION h10(n1,n2,n3,t)
!---------------------------------------------------
! hermite function basis h10
! Used for cubic Hermite spline
!---------------------------------------------------
  IMPLICIT NONE

  integer, parameter :: dp=8
  integer :: n1,n2,n3
  REAL(kind=dp), intent(in) :: t(n1,n2,n3)
  REAL(kind=dp) :: h10(n1,n2,n3)

  h10 = t**3 - 2._DP*t**2 + t
  return
END FUNCTION h10

FUNCTION h01(n1,n2,n3,t)
!---------------------------------------------------
! hermite function basis h10
! Used for cubic Hermite spline
!---------------------------------------------------
  IMPLICIT NONE

  integer, parameter :: dp=8
  integer :: n1,n2,n3
  REAL(kind=dp),intent(in) :: t(n1,n2,n3)
  REAL(kind=dp) :: h01(n1,n2,n3)

  h01 = -2._DP*t**3 + 3._DP*t**2 
  return 
END FUNCTION h01

FUNCTION h11(n1,n2,n3,t)
!---------------------------------------------------
! hermite function basis h10
! Used for cubic Hermite spline
!---------------------------------------------------
  IMPLICIT NONE

  integer, parameter :: dp=8
  integer :: n1,n2,n3
  REAL(kind=dp),intent(in) :: t(n1,n2,n3)
  REAL(kind=dp) :: h11(n1,n2,n3)

  h11 = t**3 - t**2 
  return 
END FUNCTION h11


FUNCTION stepfun(n1,n2,n3,array)
!------------------------------------------------------------------------------
! Array is an input array, 
! stepfun(i,j,k)=1 if array(i,j,k)>=0, otherwise stepfun(i,j,k)=0
!
!-------------------------------------------------------------------------------
! Created by Chen Huang (2008-Mar-16)
!--------------------------------------------------------------------------------
  IMPLICIT NONE
  
  integer, parameter :: dp=8
  integer,intent(in) :: n1,n2,n3
  integer :: ii,jj,kk
  REAL(kind=DP), intent(in) :: array(n1,n2,n3)
  REAL(kind=DP) :: stepfun(n1,n2,n3)

  stepfun = 1._DP
  do kk=1,n3
    do jj=1,n2
      do ii=1,n1
        if (array(ii,jj,kk)<0.d0) stepfun(ii,jj,kk)=0.d0
      enddo
    enddo
  enddo

  RETURN
END FUNCTION stepfun

!========================
! round the number
!========================
function  round(a)
 implicit none
 real(kind=8)  :: a
 integer :: round 
 round  =  floor(a - 0.5d0) + 1
 return
end function round

end module 

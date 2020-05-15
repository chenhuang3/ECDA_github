

!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gives the determinant of a 3x3 matrix (real)
!
! FROM OFDFT code
!------------------------------------------------------------------------------
SUBROUTINE Det3(M,DetReal) 
  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!  
  integer,parameter :: dp=8   
  REAL(kind=DP), DIMENSION(3,3), INTENT(IN) :: &
    M                      ! The matrix
  REAL(kind=DP),INTENT(out) :: &
    DetReal                ! The answer

                    !>> INTERNAL VARIABLES <<!
                     !>> INITIALIZATION <<!
                     !>> FUNCTION BODY << 

  DetReal = M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))&
               -M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1))&
               +M(1,3)*(M(2,1)*M(3,2)-M(3,1)*M(2,2))

END SUBROUTINE Det3

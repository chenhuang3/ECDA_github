
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function returns the inverse of a real 3x3 matrix. It is the user's
!   responsibilty to make sure the matrix inverted is not singular. I do not
!   test for it here. You've been warned.
!
! From OFDFT code
!------------------------------------------------------------------------------
SUBROUTINE Inverse3(M,InverseReal) 
  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<
  integer,parameter :: dp=8   
  REAL(kind=DP), DIMENSION(3,3), INTENT(IN) :: &
    M                       ! the matrix

  REAL(kind=DP), DIMENSION(3,3), INTENT(out) :: &
    InverseReal             ! the answer

                    !>> INTERNAL VARIABLES <<!  
  
  REAL(kind=DP) :: &
    d                       ! d is the determinant of M.

  INTEGER :: &
    i, j                    ! Dummy indexes

                     !>> INITIALIZATION <<!
  call Det3(M,d)
                     !>> FUNCTION BODY << 
  DO i=1, 3
    DO j=1,3
      InverseReal(i,j)=(-1)**(REAL(i,kind=DP)+REAL(j,kind=DP))/d
    END DO 
  END DO
  InverseReal(1,1)=InverseReal(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))
  InverseReal(2,1)=InverseReal(2,1)*(M(1,2)*M(3,3)-M(1,3)*M(3,2))
  InverseReal(3,1)=InverseReal(3,1)*(M(1,2)*M(2,3)-M(2,2)*M(1,3))
  InverseReal(1,2)=InverseReal(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1))
  InverseReal(2,2)=InverseReal(2,2)*(M(1,1)*M(3,3)-M(1,3)*M(3,1))
  InverseReal(3,2)=InverseReal(3,2)*(M(1,1)*M(2,3)-M(2,1)*M(1,3))
  InverseReal(1,3)=InverseReal(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
  InverseReal(2,3)=InverseReal(2,3)*(M(1,1)*M(3,2)-M(3,1)*M(1,2))
  InverseReal(3,3)=InverseReal(3,3)*(M(1,1)*M(2,2)-M(2,1)*M(1,2))
  InverseReal = TRANSPOSE(InverseReal)

END SUBROUTINE Inverse3

MODULE SWAP
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64, INT32
  IMPLICIT NONE

CONTAINS
  SUBROUTINE SWAPI64(V1, V2)
    INTEGER(KIND=INT64), INTENT(INOUT) :: V1, V2
    ! Local temp
    INTEGER(KIND=INT64) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAPI64

  SUBROUTINE SWAPI32(V1, V2)
    INTEGER(KIND=INT32), INTENT(INOUT) :: V1, V2
    ! Local temp
    INTEGER(KIND=INT32) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAPI32

  SUBROUTINE SWAPR64(V1, V2)
    REAL(KIND=REAL64), INTENT(INOUT) :: V1, V2
    ! Local temp
    REAL(KIND=REAL64) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAPR64
END MODULE SWAP

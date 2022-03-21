! Implementation of useful BINARY HEAP utilties.
MODULE BINARY_HEAP
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  IMPLICIT NONE

  INTEGER, PARAMETER :: RT=REAL64 ! <- real type
  INTEGER, PARAMETER :: IT=INT64  ! <- integer type

CONTAINS

  ! Insert a new element into a min-heap.
  SUBROUTINE INSERT_MIN(N, ARRAY, ORDER, V)
    INTEGER(KIND=IT), INTENT(INOUT) :: N
    REAL(KIND=RT),    INTENT(INOUT), DIMENSION(:) :: ARRAY
    INTEGER(KIND=IT), INTENT(INOUT), DIMENSION(:) :: ORDER
    REAL(KIND=RT),    INTENT(IN) :: V
    ! Temparary index.
    INTEGER(KIND=IT) :: I, IP
    ! Insert at the end.
    N = N + 1
    ARRAY(N) = V
    ORDER(N) = N
    ! Compare upwards to maintain HEAP property.
    I = N
    CHECK_MIN_HEAP : DO
       IP = I / 2
       IF (ARRAY(ORDER(IP)) .GT. V) THEN
          CALL SWAP_REAL(ARRAY(IP), ARRAY(I))
          CALL SWAP_INT(ORDER(IP), ORDER(I))
          I = IP
       ELSE
          EXIT CHECK_MIN_HEAP
       END IF
    END DO CHECK_MIN_HEAP
  END SUBROUTINE INSERT_MIN

  ! Pop the root from a min-heap.
  SUBROUTINE POP_MIN(N, ARRAY, ORDER, V)
    INTEGER(KIND=IT), INTENT(INOUT) :: N
    REAL(KIND=RT),    INTENT(INOUT), DIMENSION(:) :: ARRAY
    INTEGER(KIND=IT), INTENT(INOUT), DIMENSION(:) :: ORDER
    REAL(KIND=RT),    INTENT(OUT) :: V
    ! Temporary index.
    INTEGER(KIND=IT) :: I, I1, I2
    ! Get the minimum element (swap with last).
    V = ARRAY(ORDER(1))
    CALL SWAP_REAL(ARRAY(1), ARRAY(N))
    CALL SWAP_INT(ORDER(1), ORDER(N))
    N = N - 1
    ! Comparex downwards to maintain HEAP property.
    I = 1
    CHECK_MIN_HEAP : DO
       ! Get the index ofthe smaller child into I1.
       I1 = 2*I
       I2 = I1 + 1
       IF (I1 .GT. SIZE(ARRAY)) THEN
          EXIT CHECK_MIN_HEAP
       ELSE IF (I2 .GT. SIZE(ARRAY)) THEN
          I1 = I2
       ELSE IF (ARRAY(ORDER(I2)) .LT. ARRAY(ORDER(I1))) THEN
          I1 = I2
       END IF
       ! Swap down if the smaller child is less than the current.
       IF (ARRAY(ORDER(I1)) .LT. ARRAY(ORDER(I))) THEN
          CALL SWAP_REAL(ARRAY(I1), ARRAY(I))
          CALL SWAP_INT(ORDER(I1), ORDER(I))
          I = I1
       ELSE
          EXIT CHECK_MIN_HEAP
       END IF
    END DO CHECK_MIN_HEAP
  END SUBROUTINE POP_MIN

  ! Swap two integers through a temporary variable.
  SUBROUTINE SWAP_INT(V1, V2)
    INTEGER(KIND=IT), INTENT(INOUT) :: V1, V2
    INTEGER(KIND=IT) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAP_INT

  ! Swap two reals through a temporary variable.
  SUBROUTINE SWAP_REAL(V1, V2)
    REAL(KIND=RT), INTENT(INOUT) :: V1, V2
    REAL(KIND=RT) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAP_REAL

END MODULE BINARY_HEAP

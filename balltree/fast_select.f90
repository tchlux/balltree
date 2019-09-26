MODULE FAST_SELECT
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  USE SWAP
  IMPLICIT NONE

CONTAINS
  ! ------------------------------------------------------------------
  !                       FastSelect method
  ! 
  ! Given VALUES list of numbers, rearrange the elements of VALUES
  ! such that the element at index K has rank K (holds its same
  ! location as if all of VALUES were sorted). Symmetrically rearrange
  ! array INDICES to keep track of prior indices.
  ! 
  ! This algorithm uses the same conceptual approach as Floyd-Rivest,
  ! but instead of standard-deviation based selection of bounds for
  ! recursion, a rank-based method is used to pick the subset of
  ! values that is searched. This simplifies the code and improves
  ! interpretability, while achieving the same tunable performance.
  ! 
  ! Arguments:
  ! 
  !   VALUES   --  A 1D array of real numbers.
  !   INDICES  --  A 1D array of original indices for elements of VALUES.
  !   K        --  A positive integer for the rank index about which
  !                VALUES should be rearranged.
  ! Optional:
  ! 
  !   DIVISOR  --  A positive integer >= 2 that represents the
  !                division factor used for large VALUES arrays.
  !   MAX_SIZE --  An integer >= DIVISOR that represents the largest
  !                sized VALUES for which the worst-case pivot value
  !                selection is tolerable. A worst-case pivot causes
  !                O( SIZE(VALUES)^2 ) runtime. This value should be
  !                determined heuristically based on compute hardware.
  ! 
  ! Output:
  ! 
  !   The elements of the array VALUES are rearranged such that the
  !   element at position VALUES(K) is in the same location it would
  !   be if all of VALUES were in sorted order. Also known as,
  !   VALUES(K) has rank K.
  ! 
  RECURSIVE SUBROUTINE ARGSELECT_R64(VALUES, INDICES, K, DIVISOR, MAX_SIZE)
    ! Arguments
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:) :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: INDICES
    INTEGER(KIND=INT64), INTENT(IN)                  :: K
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL        :: DIVISOR, MAX_SIZE
    ! Locals
    INTEGER(KIND=INT64) :: LEFT, RIGHT, L, R, MS, D
    REAL(KIND=REAL64) :: P
    ! Initialize the divisor (for making subsets).
    IF (PRESENT(DIVISOR)) THEN ; D = DIVISOR
    ELSE IF (SIZE(VALUES) .GE. 2**23) THEN ; D = 2**5
    ELSE IF (SIZE(VALUES) .GE. 2**20) THEN ; D = 2**3
    ELSE                                   ; D = 2**2
    END IF
    ! Initialize the max size (before subsets are created).
    IF (PRESENT(MAX_SIZE)) THEN ; MS = MAX_SIZE
    ELSE                        ; MS = 2**10
    END IF
    ! Initialize LEFT and RIGHT to be the entire array.
    LEFT = 1
    RIGHT = SIZE(VALUES)
    ! Loop until done finding the K-th element.
    DO WHILE (LEFT .LT. RIGHT)
       ! Use SELECT recursively to improve the quality of the
       ! selected pivot value for large arrays.
       IF (RIGHT - LEFT .GT. MS) THEN
          ! Compute how many elements should be left and right of K
          ! to maintain the same percentile in a subset.
          L = K - K / D
          R = L + (SIZE(VALUES) / D)
          ! Perform fast select on an array a fraction of the size about K.
          CALL ARGSELECT_R64(VALUES(L:R), INDICES(L:R), K - L + 1, DIVISOR, MAX_SIZE)
       END IF
       ! Pick a partition element at position K.
       P = VALUES(K)
       L = LEFT
       R = RIGHT
       ! Move the partition element to the front of the list.
       CALL SWAPR64(VALUES(LEFT), VALUES(K))
       CALL SWAPI64(INDICES(LEFT), INDICES(K))
       ! Pre-swap the left and right elements (temporarily putting a
       ! larger element on the left) before starting the partition loop.
       IF (VALUES(RIGHT) .GT. P) THEN
          CALL SWAPR64(VALUES(LEFT), VALUES(RIGHT))
          CALL SWAPI64(INDICES(LEFT), INDICES(RIGHT))
       END IF
       ! Now partition the elements about the pivot value "T".
       DO WHILE (L .LT. R)
          CALL SWAPR64(VALUES(L), VALUES(R))
          CALL SWAPI64(INDICES(L), INDICES(R))
          L = L + 1
          R = R - 1
          DO WHILE (VALUES(L) .LT. P) ; L = L + 1 ; END DO
          DO WHILE (VALUES(R) .GT. P) ; R = R - 1 ; END DO
       END DO
       ! Place the pivot element back into its appropriate place.
       IF (VALUES(LEFT) .EQ. P) THEN
          CALL SWAPR64(VALUES(LEFT), VALUES(R))
          CALL SWAPI64(INDICES(LEFT), INDICES(R))
       ELSE
          R = R + 1
          CALL SWAPR64(VALUES(R), VALUES(RIGHT))
          CALL SWAPI64(INDICES(R), INDICES(RIGHT))
       END IF
       ! adjust left and right towards the boundaries of the subset
       ! containing the (k - left + 1)th smallest element
       IF (R .LE. K) LEFT = R + 1
       IF (K .LE. R) RIGHT = R - 1
    END DO
  END SUBROUTINE ARGSELECT_R64

  ! ------------------------------------------------------------------

  RECURSIVE SUBROUTINE ARGSELECT_I64(VALUES, INDICES, K, DIVISOR, MAX_SIZE)
    ! Arguments
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: INDICES
    INTEGER(KIND=INT64), INTENT(IN)                  :: K
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL        :: DIVISOR, MAX_SIZE
    ! Locals
    INTEGER(KIND=INT64) :: LEFT, RIGHT, L, R, MS, D, P
    ! Initialize the divisor (for making subsets).
    IF (PRESENT(DIVISOR)) THEN ; D = DIVISOR
    ELSE IF (SIZE(VALUES) .GE. 2**23) THEN ; D = 2**5
    ELSE IF (SIZE(VALUES) .GE. 2**20) THEN ; D = 2**3
    ELSE                                   ; D = 2**2
    END IF
    ! Initialize the max size (before subsets are created).
    IF (PRESENT(MAX_SIZE)) THEN ; MS = MAX_SIZE
    ELSE                        ; MS = 2**10
    END IF
    ! Initialize LEFT and RIGHT to be the entire array.
    LEFT = 1
    RIGHT = SIZE(VALUES)
    ! Loop until done finding the K-th element.
    DO WHILE (LEFT .LT. RIGHT)
       ! Use SELECT recursively to improve the quality of the
       ! selected pivot value for large arrays.
       IF (RIGHT - LEFT .GT. MS) THEN
          ! Compute how many elements should be left and right of K
          ! to maintain the same percentile in a subset.
          L = K - K / D
          R = L + (SIZE(VALUES) / D)
          ! Perform fast select on an array a fraction of the size about K.
          CALL ARGSELECT_I64(VALUES(L:R), INDICES(L:R), K - L + 1, DIVISOR, MAX_SIZE)
       END IF
       ! Pick a partition element at position K.
       P = VALUES(K)
       L = LEFT
       R = RIGHT
       ! Move the partition element to the front of the list.
       CALL SWAPI64(VALUES(LEFT), VALUES(K))
       CALL SWAPI64(INDICES(LEFT), INDICES(K))
       ! Pre-swap the left and right elements (temporarily putting a
       ! larger element on the left) before starting the partition loop.
       IF (VALUES(RIGHT) .GT. P) THEN
          CALL SWAPI64(VALUES(LEFT), VALUES(RIGHT))
          CALL SWAPI64(INDICES(LEFT), INDICES(RIGHT))
       END IF
       ! Now partition the elements about the pivot value "T".
       DO WHILE (L .LT. R)
          CALL SWAPI64(VALUES(L), VALUES(R))
          CALL SWAPI64(INDICES(L), INDICES(R))
          L = L + 1
          R = R - 1
          DO WHILE (VALUES(L) .LT. P) ; L = L + 1 ; END DO
          DO WHILE (VALUES(R) .GT. P) ; R = R - 1 ; END DO
       END DO
       ! Place the pivot element back into its appropriate place.
       IF (VALUES(LEFT) .EQ. P) THEN
          CALL SWAPI64(VALUES(LEFT), VALUES(R))
          CALL SWAPI64(INDICES(LEFT), INDICES(R))
       ELSE
          R = R + 1
          CALL SWAPI64(VALUES(R), VALUES(RIGHT))
          CALL SWAPI64(INDICES(R), INDICES(RIGHT))
       END IF
       ! adjust left and right towards the boundaries of the subset
       ! containing the (k - left + 1)th smallest element
       IF (R .LE. K) LEFT = R + 1
       IF (K .LE. R) RIGHT = R - 1
    END DO
  END SUBROUTINE ARGSELECT_I64
  ! ------------------------------------------------------------------

END MODULE FAST_SELECT

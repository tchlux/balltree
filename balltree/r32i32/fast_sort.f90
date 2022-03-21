MODULE FAST_SORT
  USE ISO_FORTRAN_ENV, ONLY: REAL32, INT32
  USE SWAP, ONLY: SWAP_I32, SWAP_R32
  IMPLICIT NONE

CONTAINS
  ! ------------------------------------------------------------------
  !                        FastSort method
  ! 
  ! This routine uses a combination of QuickSort (with modestly
  ! intelligent pivot selection) and Insertion Sort (for small arrays)
  ! to achieve very fast average case sort times for both random and
  ! partially sorted data. The pivot is selected for QuickSort as the
  ! median of the first, middle, and last values in the array.
  ! 
  ! Arguments:
  ! 
  !   VALUES   --  A 1D array of real numbers.
  !   INDICES  --  A 1D array of original indices for elements of VALUES.
  ! 
  ! Optional:
  ! 
  !   MIN_SIZE --  An positive integer that represents the largest
  !                sized VALUES for which a partition about a pivot
  !                is used to reduce the size of a an unsorted array.
  !                Any size less than this will result in the use of
  !                INSERTION_ARGSORT instead of ARGPARTITION.
  ! 
  ! Output:
  ! 
  !   The elements of the array VALUES are sorted and all elements of
  !   INDICES are sorted symmetrically (given INDICES = 1, ...,
  !   SIZE(VALUES) beforehand, final INDICES will show original index
  !   of each element of VALUES before the sort operation).
  ! 
  RECURSIVE SUBROUTINE ARGSORT(VALUES, INDICES, MIN_SIZE)
    REAL(KIND=REAL32),   INTENT(INOUT), DIMENSION(:)            :: VALUES
    INTEGER(KIND=INT32), INTENT(INOUT), DIMENSION(SIZE(VALUES)) :: INDICES
    INTEGER(KIND=INT32), INTENT(IN), OPTIONAL                   :: MIN_SIZE
    ! Local variables
    INTEGER(KIND=INT32) :: I, MS
    IF (PRESENT(MIN_SIZE)) THEN ; MS = MIN_SIZE
    ELSE                        ; MS = 64
    END IF
    ! Base case, return.
    IF (SIZE(VALUES) .LT. MS) THEN
       CALL INSERTION_ARGSORT(VALUES, INDICES)
       ! Call this function recursively after pivoting about the median.
    ELSE
       ! ---------------------------------------------------------------
       ! If you are having slow runtime with the selection of pivot values 
       ! provided by ARGPARTITION, then consider using ARGSELECT instead.
       I = ARGPARTITION(VALUES, INDICES)
       ! ---------------------------------------------------------------
       ! I = SIZE(VALUES) / 2
       ! CALL ARGSELECT(VALUES, INDICES, I)
       ! ! Requires 'USE FAST_SELECT' at top of module.
       ! ---------------------------------------------------------------
       CALL ARGSORT(VALUES(:I-1), INDICES(:I-1), MS)
       CALL ARGSORT(VALUES(I+1:), INDICES(I+1:), MS)
    END IF
  END SUBROUTINE ARGSORT

  ! This function efficiently partitions values based on the median
  ! of the first, middle, and last elements of the VALUES array. This
  ! function returns the index of the pivot.
  FUNCTION ARGPARTITION(VALUES, INDICES) RESULT(LEFT)
    REAL(KIND=REAL32),   INTENT(INOUT), DIMENSION(:)            :: VALUES
    INTEGER(KIND=INT32), INTENT(INOUT), DIMENSION(SIZE(VALUES)) :: INDICES
    INTEGER(KIND=INT32) :: LEFT, MID, RIGHT
    REAL(KIND=REAL32)   :: PIVOT
    ! Use the median of the first, middle, and last element as the
    ! pivot. Place the pivot at the end of the array.
    MID = (1 + SIZE(VALUES)) / 2
    ! Swap the first and last elements (if the last is smaller).
    IF (VALUES(SIZE(VALUES)) < VALUES(1)) THEN
       CALL SWAP_R32(VALUES(1),  VALUES(SIZE(VALUES)))
       CALL SWAP_I32(INDICES(1), INDICES(SIZE(VALUES)))
    END IF
    ! Swap the middle and first elements (if the middle is smaller).
    IF (VALUES(MID) < VALUES(SIZE(VALUES))) THEN
       CALL SWAP_R32(VALUES(MID),  VALUES(SIZE(VALUES)))
       CALL SWAP_I32(INDICES(MID), INDICES(SIZE(VALUES)))       
       ! Swap the last and first elements (if the last is smaller).
       IF (VALUES(SIZE(VALUES)) < VALUES(1)) THEN
          CALL SWAP_R32(VALUES(1),  VALUES(SIZE(VALUES)))
          CALL SWAP_I32(INDICES(1), INDICES(SIZE(VALUES)))
       END IF
    END IF
    ! Set the pivot, LEFT index and RIGHT index (skip the smallest,
    ! which is in location 1, and the pivot at the end).
    PIVOT = VALUES(SIZE(VALUES))
    LEFT  = 2
    RIGHT = SIZE(VALUES) - 1
    ! Partition all elements to the left and right side of the pivot
    ! (left if they are smaller, right if they are bigger).
    DO WHILE (LEFT < RIGHT)
       ! Loop left until we find a value that is greater or equal to pivot.
       DO WHILE (VALUES(LEFT) < PIVOT)
          LEFT = LEFT + 1
       END DO
       ! Loop right until we find a value that is less or equal to pivot (or LEFT).
       DO WHILE (RIGHT .NE. LEFT)
          IF (VALUES(RIGHT) .LT. PIVOT) EXIT
          RIGHT = RIGHT - 1
       END DO
       ! Now we know that [VALUES(RIGHT) < PIVOT < VALUES(LEFT)], so swap them.
       CALL SWAP_R32(VALUES(LEFT),  VALUES(RIGHT))
       CALL SWAP_I32(INDICES(LEFT), INDICES(RIGHT))
    END DO
    ! The last swap was done even though LEFT == RIGHT, we need to undo.
    CALL SWAP_R32(VALUES(LEFT),  VALUES(RIGHT))
    CALL SWAP_I32(INDICES(LEFT), INDICES(RIGHT))
    ! Finally, we put the pivot back into its proper location.
    CALL SWAP_R32(VALUES(LEFT),  VALUES(SIZE(VALUES)))
    CALL SWAP_I32(INDICES(LEFT), INDICES(SIZE(VALUES)))
  END FUNCTION ARGPARTITION

  ! Insertion sort (best for small lists).
  SUBROUTINE INSERTION_ARGSORT(VALUES, INDICES)
    REAL(KIND=REAL32),   INTENT(INOUT), DIMENSION(:)            :: VALUES
    INTEGER(KIND=INT32), INTENT(INOUT), DIMENSION(SIZE(VALUES)) :: INDICES
    ! Local variables.
    REAL(KIND=REAL32)   :: TEMP_VAL
    INTEGER(KIND=INT32) :: I, BEFORE, AFTER, TEMP_IND
    ! Return for the base case.
    IF (SIZE(VALUES) .LE. 1) RETURN
    ! Put the smallest value at the front of the list.
    I = MINLOC(VALUES,1)
    CALL SWAP_R32(VALUES(1),  VALUES(I))
    CALL SWAP_I32(INDICES(1), INDICES(I))
    ! Insertion sort the rest of the array.
    DO I = 3, SIZE(VALUES)
       TEMP_VAL = VALUES(I)
       TEMP_IND = INDICES(I)
       ! Search backwards in the list, 
       BEFORE = I - 1
       AFTER  = I
       DO WHILE (TEMP_VAL .LT. VALUES(BEFORE))
          VALUES(AFTER)  = VALUES(BEFORE)
          INDICES(AFTER) = INDICES(BEFORE)
          BEFORE = BEFORE - 1
          AFTER  = AFTER - 1
       END DO
       ! Put the value into its place (where it is greater than the
       ! element before it, but less than all values after it).
       VALUES(AFTER)  = TEMP_VAL
       INDICES(AFTER) = TEMP_IND
    END DO
  END SUBROUTINE INSERTION_ARGSORT
  ! ------------------------------------------------------------------

END MODULE FAST_SORT

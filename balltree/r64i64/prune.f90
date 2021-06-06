MODULE PRUNE
  USE ISO_FORTRAN_ENV, ONLY: INT64
  IMPLICIT NONE

CONTAINS

  ! Get the indices of all the inner children in this tree (the 50th
  ! percentile points) up to level 'LEVELS'.
  RECURSIVE SUBROUTINE INNER(TREE_SIZE, LEVELS, INDICES, &
       CURRENT_LEVEL, STARTING_INDEX, NEXT_INDEX)
    INTEGER(KIND=INT64), INTENT(IN)                :: TREE_SIZE, LEVELS
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL      :: CURRENT_LEVEL, STARTING_INDEX
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: NEXT_INDEX
    ! Local variables
    INTEGER(KIND=INT64) :: SI, CL, I, TMID
    ! Assign default value for 'CURRENT_LEVEL' and 'STARTING_INDEX'.
    IF (PRESENT(CURRENT_LEVEL)) THEN
       CL = CURRENT_LEVEL
       SI = STARTING_INDEX
       I = NEXT_INDEX
    ELSE
       CL = 1
       SI = 1
       I = 0
    END IF
    ! Check stopping condition.
    IF (CL+1 .GT. LEVELS) RETURN
    ! Get the middle index of the tree.
    TMID = (TREE_SIZE + 2) / 2
    ! Add the inner to the list of kept indices.
    IF ((1 .LT. TMID) .AND. (2 .LE. TREE_SIZE)) THEN
       I = I + 1
       INDICES(I) = SI + 1
    END IF
    ! Continue by recursively adding the two children trees. Subtract
    ! one from TMID to account for the current root. Skip it in indices.
    MORE_LEVELS : IF (TMID > 1) THEN
       CALL INNER(TMID-1, LEVELS, INDICES, CL+1, SI+1, I)
       ! Check to see if another child exists in the tree, if so recurse.
       HAS_OUTER : IF (TMID .LT. TREE_SIZE) THEN
          CALL INNER(TREE_SIZE-TMID, LEVELS, INDICES, CL+1, SI+TMID, I)
       END IF HAS_OUTER
    END IF MORE_LEVELS
    ! Pass the "next index" back up the recursion stack.
    IF (PRESENT(NEXT_INDEX)) NEXT_INDEX = I
  END SUBROUTINE INNER

  ! Get the indices of all the outer children in this tree (extrema
  ! points) up to level 'LEVELS'.
  RECURSIVE SUBROUTINE OUTER(TREE_SIZE, LEVELS, INDICES, &
       CURRENT_LEVEL, STARTING_INDEX, NEXT_INDEX)
    INTEGER(KIND=INT64), INTENT(IN)                :: TREE_SIZE, LEVELS
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL      :: CURRENT_LEVEL, STARTING_INDEX
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: NEXT_INDEX
    ! Local variables
    INTEGER(KIND=INT64) :: SI, CL, I, TMID
    ! Assign default value for 'CURRENT_LEVEL' and 'STARTING_INDEX'.
    IF (PRESENT(CURRENT_LEVEL)) THEN
       CL = CURRENT_LEVEL
       SI = STARTING_INDEX
       I = NEXT_INDEX
    ELSE
       CL = 1
       SI = 1
       I = 0
    END IF
    ! Check stopping condition.
    IF (CL+1 .GT. LEVELS) RETURN
    ! Get the middle index of the tree.
    TMID = (TREE_SIZE + 2) / 2
    ! Add the outer to the list of kept indices.
    IF ((1 .LT. TMID) .AND. (TMID .LT. TREE_SIZE)) THEN
       I = I + 1
       INDICES(I) = SI + TMID
    END IF
    ! Continue by recursively adding the two children trees. Subtract
    ! one from TMID to account for the current root. Skip it in indices.
    MORE_LEVELS : IF (TMID > 1) THEN
       CALL OUTER(TMID-1, LEVELS, INDICES, CL+1, SI+1, I)
       ! Check to see if another child exists in the tree, if so recurse.
       HAS_OUTER : IF (TMID .LT. TREE_SIZE) THEN
          CALL OUTER(TREE_SIZE-TMID, LEVELS, INDICES, CL+1, SI+TMID, I)
       END IF HAS_OUTER
    END IF MORE_LEVELS
    ! Pass the "next index" back up the recursion stack.
    IF (PRESENT(NEXT_INDEX)) NEXT_INDEX = I
  END SUBROUTINE OUTER

  ! Get the indices of all the top children in this tree up to level 'LEVELS'.
  RECURSIVE SUBROUTINE TOP(TREE_SIZE, LEVELS, INDICES, &
       CURRENT_LEVEL, STARTING_INDEX, NEXT_INDEX)
    INTEGER(KIND=INT64), INTENT(IN)                :: TREE_SIZE, LEVELS
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL      :: CURRENT_LEVEL, STARTING_INDEX
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: NEXT_INDEX
    ! Local variables
    INTEGER(KIND=INT64) :: SI, CL, I, TMID
    ! Assign default value for 'CURRENT_LEVEL' and 'STARTING_INDEX'.
    IF (PRESENT(CURRENT_LEVEL)) THEN
       CL = CURRENT_LEVEL
       SI = STARTING_INDEX
       I = NEXT_INDEX
    ELSE
       CL = 1
       SI = 1
       I = 0
    END IF
    ! Check stopping condition.
    IF (CL .GT. LEVELS) RETURN
    ! Get the middle index of the tree.
    TMID = (TREE_SIZE + 2) / 2
    ! Add the root to the list of kept indices.
    I = I + 1
    INDICES(I) = SI
    ! Continue by recursively adding the two children trees. Subtract
    ! one from TMID to account for the current root. Skip it in indices.
    MORE_LEVELS : IF (TMID > 1) THEN
       CALL TOP(TMID-1, LEVELS, INDICES, CL+1, SI+1, I)
       ! Check to see if another child exists in the tree, if so recurse.
       HAS_OUTER : IF (TMID .LT. TREE_SIZE) THEN
          CALL TOP(TREE_SIZE-TMID, LEVELS, INDICES, CL+1, SI+TMID, I)
       END IF HAS_OUTER
    END IF MORE_LEVELS
    ! Pass the "next index" back up the recursion stack.
    IF (PRESENT(NEXT_INDEX)) NEXT_INDEX = I
  END SUBROUTINE TOP

  ! Get the indices of the level "LAYER" of this tree, with 0 being
  ! the root. List of indices is returned in "INDICES".
  RECURSIVE SUBROUTINE LEVEL(TREE_SIZE, LAYER, INDICES, &
       FOUND, STARTING_INDEX)
    INTEGER(KIND=INT64), INTENT(IN)                :: TREE_SIZE, LAYER
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    INTEGER(KIND=INT64), INTENT(INOUT)             :: FOUND
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL :: STARTING_INDEX
    ! Local variables
    INTEGER(KIND=INT64) :: SI, TMID
    IF (PRESENT(STARTING_INDEX)) THEN
       SI = STARTING_INDEX
    ELSE
       SI = 1
    END IF
    ! Check stopping condition.
    IF (LAYER .EQ. 0) THEN
       ! Add this node to the list (if it is the right level.
       FOUND = FOUND + 1
       INDICES(FOUND) = SI
    ELSE
       ! Get the middle index of the tree.
       TMID = (TREE_SIZE + 2) / 2
       ! Continue by recursively adding the two children trees. Subtract
       ! one from TMID to account for the current root. Skip it in indices.
       MORE_LEVELS : IF (TMID > 1) THEN
          CALL LEVEL(TMID-1, LAYER-1, INDICES, FOUND, SI+1)
          ! Check to see if another child exists in the tree, if so recurse.
          HAS_OUTER : IF (TMID .LT. TREE_SIZE) THEN
             CALL LEVEL(TREE_SIZE-TMID, LAYER-1, INDICES, FOUND, SI+TMID)
          END IF HAS_OUTER
       END IF MORE_LEVELS
    END IF
  END SUBROUTINE LEVEL


END MODULE PRUNE

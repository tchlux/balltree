MODULE BALL_TREE_R64
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  USE SWAP,        ONLY: SWAP_I64
  USE FAST_SELECT, ONLY: ARGSELECT_R64
  USE FAST_SORT,   ONLY: ARGSORT
  IMPLICIT NONE

CONTAINS

  ! Use the 1-norm to compute the distance between two real vectors.
  FUNCTION DISTANCE(X, Y) RESULT(D)
    REAL(KIND=REAL64), DIMENSION(:), INTENT(IN) :: X, Y
    REAL(KIND=REAL64) :: D
    ! D = SUM(ABS(X(:) - Y(:)))
    D = SQRT(SUM(X(:)**2) + SUM(Y(:)**2) - 2 * DOT_PRODUCT(X(:),Y(:)))
  END FUNCTION DISTANCE

  ! Re-arrange elements of POINTS into a binary ball tree.
  RECURSIVE SUBROUTINE BUILD_TREE(POINTS, RADII, ORDER, ROOT, LEAF_SIZE)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:,:) :: POINTS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:) :: RADII
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL :: ROOT, LEAF_SIZE
    ! Local variables
    INTEGER(KIND=INT64) :: CENTER_IDX, MID, I, J, LS
    REAL(KIND=REAL64), DIMENSION(SIZE(POINTS,1)) :: PT
    REAL(KIND=REAL64), DIMENSION(SIZE(ORDER))    :: DISTS
    REAL(KIND=REAL64) :: MAX_DIST, DIST
    ! Set the leaf size to 1 by default (most possible work required,
    ! but guarantees successful use with any leaf size).
    IF (PRESENT(LEAF_SIZE)) THEN ; LS = LEAF_SIZE
    ELSE                         ; LS = 1
    END IF
    ! Set the index of the 'split' 
    IF (PRESENT(ROOT)) THEN ; CENTER_IDX = ROOT
    ELSE
       ! 1) Compute distances between first point (random) and all others.
       ! 2) Pick the furthest point (on conv hull) from first as the center node.
       J = ORDER(1)
       PT(:) = POINTS(:,J)
       CENTER_IDX = 1
       MAX_DIST = 0.0_REAL64
       DO I = 2, SIZE(ORDER)
          DIST = DISTANCE(POINTS(:,ORDER(I)), PT(:))
          ! If this is a new max distance, record this point.
          IF (DIST .GT. MAX_DIST) THEN
             MAX_DIST = DIST
             CENTER_IDX = I
          END IF
       END DO
       ! Now CENTER_IDX is the selected center for this node in tree.
    END IF
    ! Move the "center" to the first position.
    CALL SWAP_I64(ORDER(1), ORDER(CENTER_IDX))
    ! Measure distance beween "center" node and all other points.
    J = ORDER(1)
    PT(:) = POINTS(:,J)
    DISTS(1) = 0.0_REAL64
    CENTER_TO_ALL : DO I = 2, SIZE(ORDER)
       DISTS(I) = DISTANCE(POINTS(:,ORDER(I)), PT(:))
    END DO CENTER_TO_ALL
    ! Base case for recursion, once we have few enough points, exit.
    IF (SIZE(ORDER) .LE. LS) THEN
       RADII(ORDER(1)) = MAXVAL(DISTS(:))
       IF (SIZE(ORDER) .GT. 1) RADII(ORDER(2:)) = 0.0_REAL64
       RETURN
    ELSE IF (SIZE(ORDER) .EQ. 2) THEN
       ! If the leaf size is 1 and there are only 2 elements, store
       ! the radius and exit (since there are no further steps.
       RADII(ORDER(1)) = MAXVAL(DISTS)
       RADII(ORDER(2)) = 0.0_REAL64
       RETURN
    END IF
    ! Rearrange "DISTS" about the median value.
    ! Compute the last index that will belong "inside" this node.
    MID = (SIZE(ORDER) + 2) / 2
    CALL ARGSELECT_R64(DISTS(2:), ORDER(2:), MID - 1)
    ! Now ORDER has been rearranged such that the median distance
    ! element of POINTS is at the median location.
    ! Identify the furthest point (must be in second half of list).
    I = MID + MAXLOC(DISTS(MID+1:),1)
    ! Store the "radius" of this ball, the furthest point.
    RADII(ORDER(1)) = DISTS(I)
    ! Move the furthest point into the spot after the median (outer root).
    CALL SWAP_I64(ORDER(I), ORDER(MID+1))
    ! Move the median point (furthest "interior") to the front (inner root).
    CALL SWAP_I64(ORDER(2), ORDER(MID))
    !$OMP PARALLEL NUM_THREADS(2)
    !$OMP SECTIONS
    !$OMP SECTION
    ! Recurisively create this tree.
    !   build a tree with the root being the furthest from this center
    !   for the remaining "interior" points of this center node.
    CALL BUILD_TREE(POINTS, RADII, ORDER(2:MID), 1_INT64, LS)
    !$OMP SECTION
    !   build a tree with the root being the furthest from this center
    !   for the remaining "exterior" points of this center node.
    !   Only perform this operation if there are >0 points available.
    IF (MID < SIZE(ORDER)) &
         CALL BUILD_TREE(POINTS, RADII, ORDER(MID+1:), 1_INT64, LS)
    !$OMP END SECTIONS
    !$OMP END PARALLEL
  END SUBROUTINE BUILD_TREE

  ! Compute the K nearest elements of TREE to each point in POINTS.
  SUBROUTINE NEAREST(POINTS, K, TREE, RADII, ORDER, LEAF_SIZE, &
       INDICES, DISTS, TO_SEARCH)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: POINTS, TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: LEAF_SIZE, K
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: DISTS
    INTEGER(KIND=INT64), INTENT(IN),  OPTIONAL    :: TO_SEARCH
    ! Local variables.
    INTEGER(KIND=INT64) :: I, B, BUDGET
    INTEGER(KIND=INT64), DIMENSION(K+LEAF_SIZE+2) :: INDS_BUFFER
    REAL(KIND=REAL64),   DIMENSION(K+LEAF_SIZE+2) :: DISTS_BUFFER
    IF (PRESENT(TO_SEARCH)) THEN ; BUDGET = MAX(K, TO_SEARCH)
    ELSE ; BUDGET = SIZE(ORDER) ; END IF
    ! For each point in this set, use the recursive branching
    ! algorithm to identify the nearest elements of TREE.
    !$OMP PARALLEL DO PRIVATE(INDS_BUFFER, DISTS_BUFFER, B)
    DO I = 1, SIZE(POINTS,2)
       B = BUDGET
       CALL PT_NEAREST(POINTS(:,I), K, TREE, RADII, ORDER, &
            LEAF_SIZE, INDS_BUFFER, DISTS_BUFFER, B)
       ! Sort the first K elements of the temporary arry for return.
       INDICES(:,I) = INDS_BUFFER(:K)
       DISTS(:,I) = DISTS_BUFFER(:K)
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE NEAREST

  ! Compute the K nearest elements of TREE to each point in POINTS.
  RECURSIVE SUBROUTINE PT_NEAREST(POINT, K, TREE, RADII, ORDER, &
       LEAF_SIZE, INDICES, DISTS, CHECKS, FOUND)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: POINT
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(:) :: DISTS
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: CHECKS
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: FOUND
    ! Local variables
    INTEGER(KIND=INT64) :: F, I, I1, I2, MID, ALLOWED_CHECKS
    REAL(KIND=REAL64)   :: D, D1, D2
    ! Initialize FOUND for first call, if FOUND is present then
    ! this must not be first and all are present.
    INITIALIZE : IF (PRESENT(FOUND)) THEN
       ALLOWED_CHECKS = CHECKS
       IF (ALLOWED_CHECKS .LE. 0) RETURN
       F = FOUND
    ELSE
       ! Initialize the remaining checks to search.
       IF (PRESENT(CHECKS)) THEN ; ALLOWED_CHECKS = CHECKS - 1
       ELSE ; ALLOWED_CHECKS = SIZE(ORDER) - 1 ; END IF
       ! Start at index 0 (added onto current index). Measure distance to root.
       INDICES(1) = ORDER(1)
       DISTS(1) = DISTANCE(POINT(:), TREE(:,ORDER(1)))
       ! Set the "points found" to be 1.
       F = 1
    END IF INITIALIZE
    ! If this is NOT a leaf node, then recurse.
    BRANCH_OR_LEAF : IF (SIZE(ORDER) .GT. LEAF_SIZE) THEN
       ALLOWED_CHECKS = ALLOWED_CHECKS - 1
       ! Measure distance to inner child.
       I1 = ORDER(2)
       D1 = DISTANCE(POINT(:),TREE(:,I1))
       ! Store this distance calculation and index.
       F = F + 1
       INDICES(F) = I1
       DISTS(F) = D1
       ! Measure distance to outer child the same as above, after
       ! checking to see if there *is* an outer child.
       MID = (SIZE(ORDER) + 2) / 2
       IF (MID+1 > SIZE(ORDER)) THEN
          I2 = 0
          D2 = HUGE(D2)
       ELSE
          I2 = ORDER(MID+1)
          IF (I2 .NE. I1) THEN
             ALLOWED_CHECKS = ALLOWED_CHECKS - 1
             D2 = DISTANCE(POINT(:),TREE(:,I2))
             ! Store this distance calculation and index.
             F = F + 1
             INDICES(F) = I2
             DISTS(F) = D2
          ELSE ; D2 = HUGE(D2)
          END IF
       END IF
       ! Re-organize the list of closest points, pushing them to first K spots.
       CALL ARGSELECT_R64(DISTS(:F), INDICES(:F), MIN(K,F))
       F = MIN(K,F)
       ! Store the maximum distance.
       D = MAXVAL(DISTS(:F),1)
       ! Determine which child to search (depth-first search) based
       ! on which child region the point lands in from the root.
       INNER_CHILD_CLOSER : IF (D1 .LE. D2) THEN
          ! Search the inner child if it could contain a nearer point.
          SEARCH_INNER1 : IF ((F .LT. K) .OR. (D .GT. D1 - RADII(I1))) THEN
             CALL PT_NEAREST(POINT, K, TREE, RADII, ORDER(2:MID), &
                  LEAF_SIZE, INDICES, DISTS, ALLOWED_CHECKS, F)
          END IF SEARCH_INNER1
          ! Search the outer child if it could contain a nearer point.
          SEARCH_OUTER1 : IF (((F .LT. K) .OR. (D .GT. D2 - RADII(I2))) &
               .AND. (I2 .NE. I1) .AND. (I2 .GT. 0)) THEN
             CALL PT_NEAREST(POINT, K, TREE, RADII, ORDER(MID+1:), &
                  LEAF_SIZE, INDICES, DISTS, ALLOWED_CHECKS, F)
          END IF SEARCH_OUTER1
       ELSE
          ! Search the outer child if it could contain a nearer point.
          SEARCH_OUTER2 : IF (((F .LT. K) .OR. (D .GT. D2 - RADII(I2))) &
               .AND. (I2 .NE. I1) .AND. (I2 .GT. 0)) THEN
             CALL PT_NEAREST(POINT, K, TREE, RADII, ORDER(MID+1:), &
                  LEAF_SIZE, INDICES, DISTS, ALLOWED_CHECKS, F)
          END IF SEARCH_OUTER2
          ! Search the inner child if it could contain a nearer point.
          SEARCH_INNER2 : IF ((F .LT. K) .OR. (D .GT. D1 - RADII(I1))) THEN
             CALL PT_NEAREST(POINT, K, TREE, RADII, ORDER(2:MID), &
                  LEAF_SIZE, INDICES, DISTS, ALLOWED_CHECKS, F)
          END IF SEARCH_INNER2
       END IF INNER_CHILD_CLOSER
    ! Since this is a leaf node, we measure distance to all children.
    ELSE
       DIST_TO_CHILDREN : DO I = 2, SIZE(ORDER)
          IF (ALLOWED_CHECKS .LE. 0) EXIT DIST_TO_CHILDREN
          ALLOWED_CHECKS = ALLOWED_CHECKS - 1
          ! Measure distance to all children of this node.
          D = DISTANCE(POINT(:), TREE(:,ORDER(I)))
          ! Store this distance.
          F = F + 1
          DISTS(F) = D
          INDICES(F) = ORDER(I)
       END DO DIST_TO_CHILDREN
       ! Reduce the kept points to only those that are closest.
       CALL ARGSELECT_R64(DISTS(:F), INDICES(:F), MIN(K,F))
       F = MIN(K, F)
    END IF BRANCH_OR_LEAF

    ! Handle closing operations..
    SORT_K : IF (PRESENT(FOUND)) THEN
       ! This is not the root, we need to pass the updated value of
       ! FOUND back up the recrusion stack.
       FOUND = F
       CHECKS = ALLOWED_CHECKS
    ELSE
       ! This is the root, initial caller. Sort the distances for return.
       CALL ARGSELECT_R64(DISTS(:F), INDICES(:F), K)
       CALL ARGSORT(DISTS(:K), INDICES(:K))
    END IF SORT_K
  END SUBROUTINE PT_NEAREST

  ! Re-organize a built tree so that it is more usefully packed in memory.
  SUBROUTINE FIX_ORDER(POINTS, RADII, ORDER, COPY)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:,:) :: POINTS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:)   :: ORDER
    LOGICAL, INTENT(IN), OPTIONAL :: COPY
    LOGICAL :: SHOULD_COPY
    INTEGER(KIND=INT64) :: I
    ! Default to copy (in memory) if there is less than 1 GB of data.
    IF (PRESENT(COPY)) THEN ; SHOULD_COPY = COPY
    ELSE ;                    SHOULD_COPY = SIZEOF(POINTS) .LE. 2**30
    END IF
    ! Reorder all of the data.

    ! IF (SHOULD_COPY) THEN
    POINTS(:,:) = POINTS(:,ORDER)
    ! ELSE
    !    OPEN( 1, STATUS='SCRATCH', ACCESS='WRITE')
    !    DO I = 1, SIZE(POINTS,2)
    !    END DO
    !    ! Write the data in the correct order to a file.
    !    ! Read the correctly ordered data back into memory
    ! END IF

    ! Always copy the square sums and the radii in memory (shouldn't be bad).
    RADII(:) = RADII(ORDER)
    FORALL (I=1:SIZE(ORDER)) ORDER(I) = I
  END SUBROUTINE FIX_ORDER

END MODULE BALL_TREE_R64

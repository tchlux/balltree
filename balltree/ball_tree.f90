MODULE BALL_TREE
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64, INT8
  USE SWAP
  USE FAST_SELECT
  USE FAST_SORT
  IMPLICIT NONE

CONTAINS

  ! Re-arrange elements of POINTS into a binary ball tree.
  RECURSIVE SUBROUTINE BUILD_TREE_R64(POINTS, SQ_SUMS, RADII, ORDER,&
       ROOT, LEAF_SIZE, COMPUTED_SQ_SUMS)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:,:) :: POINTS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:) :: SQ_SUMS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:) :: RADII
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL :: ROOT, LEAF_SIZE
    LOGICAL,             INTENT(IN), OPTIONAL :: COMPUTED_SQ_SUMS
    ! Local variables
    INTEGER(KIND=INT64) :: CENTER_IDX, MID, I, J, LS
    REAL(KIND=REAL64), DIMENSION(SIZE(POINTS,1)) :: PT
    REAL(KIND=REAL64), DIMENSION(SIZE(ORDER)) :: SQ_DISTS
    REAL(KIND=REAL64) :: MAX_SQ_DIST, SQ_DIST
    ! Set the leaf size to 1 by default (most possible work required,
    ! but guarantees successful use with any leaf size).
    IF (PRESENT(LEAF_SIZE)) THEN ; LS = LEAF_SIZE
    ELSE                         ; LS = 1
    END IF
    ! If no squared sums were provided, compute them.
    IF (.NOT. PRESENT(COMPUTED_SQ_SUMS) .OR. &
         .NOT. COMPUTED_SQ_SUMS) THEN
       DO I = 1, SIZE(POINTS,2)
          SQ_SUMS(I) = SUM(POINTS(:,I)**2)
       END DO
    END IF
    ! Set the index of the 'split' 
    IF (PRESENT(ROOT)) THEN ; CENTER_IDX = ROOT
    ELSE
       ! 1) Compute distances between first point (random) and all others.
       ! 2) Pick the furthest point (on conv hull) from first as the center node.
       J = ORDER(1)
       PT(:) = POINTS(:,J)
       CENTER_IDX = 1
       MAX_SQ_DIST = 0.0_REAL64
       DO I = 2, SIZE(ORDER)
          SQ_DIST = SQ_SUMS(J) + SQ_SUMS(ORDER(I)) - &
               2 * DOT_PRODUCT(POINTS(:,ORDER(I)), PT(:))
          ! If this is a new max distance, record this point.
          IF (SQ_DIST .GT. MAX_SQ_DIST) THEN
             MAX_SQ_DIST = SQ_DIST
             CENTER_IDX = I
          END IF
       END DO
       ! Now CENTER_IDX is the selected center for this node in tree.
    END IF
    ! Move the "center" to the first position.
    CALL SWAPI64(ORDER(1), ORDER(CENTER_IDX))
    ! Measure squared distance beween "center" node and all other points.
    J = ORDER(1)
    PT(:) = POINTS(:,J)
    SQ_DISTS(1) = 0.0_REAL64
    CENTER_TO_ALL : DO I = 2, SIZE(ORDER)
       SQ_DISTS(I) = SQ_SUMS(J) + SQ_SUMS(ORDER(I)) - &
            2 * DOT_PRODUCT(POINTS(:,ORDER(I)), PT(:))
    END DO CENTER_TO_ALL
    ! Base case for recursion, once we have few enough points, exit.
    IF (SIZE(ORDER) .LE. LS) THEN
       RADII(ORDER(1)) = SQRT(MAXVAL(SQ_DISTS))
       IF (SIZE(ORDER) .GT. 1) RADII(ORDER(2:)) = 0.0_REAL64
       RETURN
    ELSE IF (SIZE(ORDER) .EQ. 2) THEN
       ! If the leaf size is 1 and there are only 2 elements, store
       ! the radius and exit (since there are no further steps.
       RADII(ORDER(1)) = SQRT(SQ_DISTS(2))
       RADII(ORDER(2)) = 0.0_REAL64
       RETURN
    END IF
    ! Rearrange "SQ_DISTS" about the median value.
    ! Compute the last index that will belong "inside" this node.
    MID = (SIZE(ORDER) + 2) / 2
    CALL ARGSELECT_R64(SQ_DISTS(2:), ORDER(2:), MID - 1)
    ! Now ORDER has been rearranged such that the median distance
    ! element of POINTS is at the median location.
    ! Identify the furthest point (must be in second half of list).
    I = MID + MAXLOC(SQ_DISTS(MID+1:),1)
    ! Store the "radius" of this ball, the furthest point.
    RADII(ORDER(1)) = SQRT(SQ_DISTS(I))
    ! Move the furthest point into the spot after the median (outer root).
    CALL SWAPI64(ORDER(I), ORDER(MID+1))
    ! Move the median point (furthest "interior") to the front (inner root).
    CALL SWAPI64(ORDER(2), ORDER(MID))
    ! Recurisively create this tree.
    !   build a tree with the root being the furthest from this center
    !   for the remaining "interior" points of this center node.
    CALL BUILD_TREE_R64(POINTS, SQ_SUMS, RADII, ORDER(2:MID), 1_INT64, LS, .TRUE.)
    !   build a tree with the root being the furthest from this center
    !   for the remaining "exterior" points of this center node.
    !   Only perform this operation if there are >0 points available.
    IF (MID < SIZE(ORDER)) &
         CALL BUILD_TREE_R64(POINTS, SQ_SUMS, RADII, &
         ORDER(MID+1:), 1_INT64, LS, .TRUE.)
  END SUBROUTINE BUILD_TREE_R64


  ! Re-arrange elements of POINTS into a binary ball tree.
  RECURSIVE SUBROUTINE BUILD_TREE_I8(POINTS, SQ_SUMS, RADII, ORDER,&
       ROOT, LEAF_SIZE, COMPUTED_SQ_SUMS)
    INTEGER(KIND=INT8),  INTENT(INOUT), DIMENSION(:,:) :: POINTS
    INTEGER(KIND=INT64), INTENT(OUT),   DIMENSION(:) :: SQ_SUMS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:) :: RADII
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL :: ROOT, LEAF_SIZE
    LOGICAL,             INTENT(IN), OPTIONAL :: COMPUTED_SQ_SUMS
    ! Local variables
    INTEGER(KIND=INT64) :: CENTER_IDX, MID, I, J, LS
    INTEGER(KIND=INT64), DIMENSION(SIZE(POINTS,1)) :: PT
    INTEGER(KIND=INT64), DIMENSION(SIZE(ORDER))    :: SQ_DISTS
    INTEGER(KIND=INT64) :: MAX_SQ_DIST, SQ_DIST
    ! Set the leaf size to 1 by default (most possible work required,
    ! but guarantees successful use with any leaf size).
    IF (PRESENT(LEAF_SIZE)) THEN ; LS = LEAF_SIZE
    ELSE                         ; LS = 1
    END IF
    ! If no squared sums were provided, compute them.
    IF (.NOT. PRESENT(COMPUTED_SQ_SUMS) .OR. &
         .NOT. COMPUTED_SQ_SUMS) THEN
       DO I = 1, SIZE(POINTS,2)
          SQ_SUMS(I) = SUM(POINTS(:,I)**2)
       END DO
    END IF
    ! Set the index of the 'split' 
    IF (PRESENT(ROOT)) THEN ; CENTER_IDX = ROOT
    ELSE
       ! 1) Compute distances between first point (random) and all others.
       ! 2) Pick the furthest point (on conv hull) from first as the center node.
       J = ORDER(1)
       PT(:) = INT(POINTS(:,J), INT64)
       CENTER_IDX = 1
       MAX_SQ_DIST = 0
       DO I = 2, SIZE(ORDER)
          SQ_DIST = SQ_SUMS(J) + SQ_SUMS(ORDER(I)) - &
               2 * DOT_PRODUCT(POINTS(:,ORDER(I)), PT(:))
          ! If this is a new max distance, record this point.
          IF (SQ_DIST .GT. MAX_SQ_DIST) THEN
             MAX_SQ_DIST = SQ_DIST
             CENTER_IDX = I
          END IF
       END DO
       ! Now CENTER_IDX is the selected center for this node in tree.
    END IF
    ! Move the "center" to the first position.
    CALL SWAPI64(ORDER(1), ORDER(CENTER_IDX))
    ! Measure squared distance beween "center" node and all other points.
    J = ORDER(1)
    PT(:) = INT(POINTS(:,J), INT64)
    SQ_DISTS(1) = 0
    CENTER_TO_ALL : DO I = 2, SIZE(ORDER)
       SQ_DISTS(I) = SQ_SUMS(J) + SQ_SUMS(ORDER(I)) - &
            2 * DOT_PRODUCT(POINTS(:,ORDER(I)), PT(:))
    END DO CENTER_TO_ALL
    ! Base case for recursion, once we have few enough points, exit.
    IF (SIZE(ORDER) .LE. LS) THEN
       RADII(ORDER(1)) = SQRT( REAL(MAXVAL(SQ_DISTS),REAL64) )
       IF (SIZE(ORDER) .GT. 1) RADII(ORDER(2:)) = 0.0_REAL64
       RETURN
    ELSE IF (SIZE(ORDER) .EQ. 2) THEN
       ! If the leaf size is 1 and there are only 2 elements, store
       ! the radius and exit (since there are no further steps.
       RADII(ORDER(1)) = SQRT( REAL(SQ_DISTS(2),REAL64) )
       RADII(ORDER(2)) = 0.0_REAL64
       RETURN
    END IF
    ! Rearrange "SQ_DISTS" about the median value.
    ! Compute the last index that will belong "inside" this node.
    MID = (SIZE(ORDER) + 2) / 2
    CALL ARGSELECT_I64(SQ_DISTS(2:), ORDER(2:), MID - 1)
    ! Now ORDER has been rearranged such that the median distance
    ! element of POINTS is at the median location.
    ! Identify the furthest point (must be in second half of list).
    I = MID + MAXLOC(SQ_DISTS(MID+1:),1)
    ! Store the "radius" of this ball, the furthest point.
    RADII(ORDER(1)) = SQRT( REAL(SQ_DISTS(I),REAL64) )
    ! Move the furthest point into the spot after the median (outer root).
    CALL SWAPI64(ORDER(I), ORDER(MID+1))
    ! Move the median point (furthest "interior") to the front (inner root).
    CALL SWAPI64(ORDER(2), ORDER(MID))
    ! Recurisively create this tree.
    !   build a tree with the root being the furthest from this center
    !   for the remaining "interior" points of this center node.
    CALL BUILD_TREE_I8(POINTS, SQ_SUMS, RADII, ORDER(2:MID), 1_INT64, LS, .TRUE.)
    !   build a tree with the root being the furthest from this center
    !   for the remaining "exterior" points of this center node.
    !   Only perform this operation if there are >0 points available.
    IF (MID < SIZE(ORDER)) &
         CALL BUILD_TREE_I8(POINTS, SQ_SUMS, RADII, &
         ORDER(MID+1:), 1_INT64, LS, .TRUE.)
  END SUBROUTINE BUILD_TREE_I8

  ! Compute the K nearest elements of TREE to each point in POINTS.
  SUBROUTINE NEAREST_R64(POINTS, K, TREE, SQ_SUMS, RADII, ORDER, &
       LEAF_SIZE, INDICES, DISTS)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: POINTS, TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: DISTS
    ! Local variables.
    INTEGER :: I
    INTEGER(KIND=INT64), DIMENSION(K+LEAF_SIZE+2) :: INDS_BUFFER
    REAL(KIND=REAL64),   DIMENSION(K+LEAF_SIZE+2) :: DISTS_BUFFER
    ! For each point in this set, use the recursive branching
    ! algorithm to identify the nearest elements of TREE.
    DO I = 1, SIZE(POINTS,2)
       CALL PT_NEAREST_R64(POINTS(:,I), K, TREE, SQ_SUMS, RADII, ORDER, &
            LEAF_SIZE, INDS_BUFFER, DISTS_BUFFER)
       ! Sort the first K elements of the temporary arry for return.
       INDICES(:,I) = INDS_BUFFER(:K)
       DISTS(:,I) = DISTS_BUFFER(:K)
    END DO
  END SUBROUTINE NEAREST_R64

  ! Compute the K nearest elements of TREE to each point in POINTS.
  SUBROUTINE NEAREST_I8(POINTS, K, TREE, SQ_SUMS, RADII, ORDER, &
       LEAF_SIZE, INDICES, DISTS)
    INTEGER(KIND=INT8),  INTENT(IN), DIMENSION(:,:) :: POINTS, TREE
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64),   INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:)   :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)                 :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: DISTS
    ! Local variables.
    INTEGER :: I
    INTEGER(KIND=INT64), DIMENSION(K+LEAF_SIZE+2) :: INDS_BUFFER
    REAL(KIND=REAL64),   DIMENSION(K+LEAF_SIZE+2) :: DISTS_BUFFER
    ! For each point in this set, use the recursive branching
    ! algorithm to identify the nearest elements of TREE.
    DO I = 1, SIZE(POINTS,2)
       CALL PT_NEAREST_I8(INT(POINTS(:,I),INT64), K, TREE, SQ_SUMS, &
            RADII, ORDER, LEAF_SIZE, INDS_BUFFER, DISTS_BUFFER)
       ! Sort the first K elements of the temporary arry for return.
       INDICES(:,I) = INDS_BUFFER(:K)
       DISTS(:,I) = DISTS_BUFFER(:K)
    END DO
  END SUBROUTINE NEAREST_I8

  ! Compute the K nearest elements of TREE to each point in POINTS.
  RECURSIVE SUBROUTINE PT_NEAREST_R64(POINT, K, TREE, SQ_SUMS, RADII, &
       ORDER, LEAF_SIZE, INDICES, DISTS, FOUND, PT_SS)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: POINT
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(:) :: DISTS
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: FOUND
    REAL(KIND=REAL64),   INTENT(IN),    OPTIONAL   :: PT_SS
    ! Local variables
    INTEGER(KIND=INT64) :: F, I, I1, I2, MID
    REAL(KIND=REAL64)   :: D, D1, D2
    REAL(KIND=REAL64)   :: PS
    ! Initialize FOUND for first call, if FOUND is present then
    ! this must not be first and all are present.
    INITIALIZE : IF (PRESENT(FOUND)) THEN
       F = FOUND
       PS = PT_SS
    ELSE
       ! Start at index 0 (added onto current index). Compute squared sum.
       PS = SUM(POINT(:)**2)
       ! Measure distance to root.
       INDICES(1) = ORDER(1)
       DISTS(1) = SQRT(PS + SQ_SUMS(ORDER(1)) - &
            2*DOT_PRODUCT(POINT(:), TREE(:,ORDER(1))))
       ! Set the "points found" to be 1.
       F = 1
    END IF INITIALIZE
    ! If this is NOT a leaf node, then recurse.
    BRANCH_OR_LEAF : IF (SIZE(ORDER) .GT. LEAF_SIZE) THEN
       ! Measure distance to inner child.
       I1 = ORDER(2)
       D1 = SQRT(PS + SQ_SUMS(I1) - &
            2*DOT_PRODUCT(POINT(:),TREE(:,I1)))
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
             D2 = SQRT(PS + SQ_SUMS(I2) - &
                  2*DOT_PRODUCT(POINT(:),TREE(:,I2)))
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
             CALL PT_NEAREST_R64(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(2:MID), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_INNER1
          ! Search the outer child if it could contain a nearer point.
          SEARCH_OUTER1 : IF (((F .LT. K) .OR. (D .GT. D2 - RADII(I2))) &
               .AND. (I2 .NE. I1) .AND. (I2 .GT. 0)) THEN
             CALL PT_NEAREST_R64(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(MID+1:), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_OUTER1
       ELSE
          ! Search the outer child if it could contain a nearer point.
          SEARCH_OUTER2 : IF (((F .LT. K) .OR. (D .GT. D2 - RADII(I2))) &
               .AND. (I2 .NE. I1) .AND. (I2 .GT. 0)) THEN
             CALL PT_NEAREST_R64(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(MID+1:), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_OUTER2
          ! Search the inner child if it could contain a nearer point.
          SEARCH_INNER2 : IF ((F .LT. K) .OR. (D .GT. D1 - RADII(I1))) THEN
             CALL PT_NEAREST_R64(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(2:MID), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_INNER2
       END IF INNER_CHILD_CLOSER

    ! Since this is a leaf node, we measure distance to all children.
    ELSE
       DIST_TO_CHILDREN : DO I = 2, SIZE(ORDER)
          ! Measure distance to all children of this node.
          D = SQRT(PS + SQ_SUMS(ORDER(I)) - &
               2*DOT_PRODUCT(POINT(:),TREE(:,ORDER(I))))
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
    ELSE
       ! This is the root, initial caller. Sort the distances for return.
       CALL ARGSELECT_R64(DISTS(:F), INDICES(:F), K)
       CALL ARGSORT(DISTS(:K), INDICES(:K))
    END IF SORT_K
  END SUBROUTINE PT_NEAREST_R64


  ! Compute the K nearest elements of TREE to each point in POINTS.
  RECURSIVE SUBROUTINE PT_NEAREST_I8(POINT, K, TREE, SQ_SUMS, RADII, &
       ORDER, LEAF_SIZE, INDICES, DISTS, FOUND, PT_SS)
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:)   :: POINT
    INTEGER(KIND=INT8),  INTENT(IN), DIMENSION(:,:) :: TREE
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64),   INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:)   :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)                 :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(:) :: DISTS
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: FOUND
    INTEGER(KIND=INT64), INTENT(IN),    OPTIONAL   :: PT_SS
    ! Local variables
    INTEGER(KIND=INT64) :: F, I, I1, I2, MID
    REAL(KIND=REAL64)   :: D, D1, D2
    INTEGER(KIND=INT64) :: PS
    ! Initialize FOUND for first call, if FOUND is present then
    ! this must not be first and all are present.
    INITIALIZE : IF (PRESENT(FOUND)) THEN
       F = FOUND
       PS = PT_SS
    ELSE
       ! Start at index 0 (added onto current index). Compute squared sum.
       PS = SUM(POINT(:)**2)
       ! Measure distance to root.
       INDICES(1) = ORDER(1)
       DISTS(1) = SQRT( REAL(PS + SQ_SUMS(ORDER(1)) - &
            2*DOT_PRODUCT(POINT(:), TREE(:,ORDER(1))), REAL64) )
       ! Set the "points found" to be 1.
       F = 1
    END IF INITIALIZE
    ! If this is NOT a leaf node, then recurse.
    BRANCH_OR_LEAF : IF (SIZE(ORDER) .GT. LEAF_SIZE) THEN
       ! Measure distance to inner child.
       I1 = ORDER(2)
       D1 = SQRT( REAL(PS + SQ_SUMS(I1) - &
            2*DOT_PRODUCT(POINT(:),TREE(:,I1)), REAL64) )
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
          I2 = ORDER(MID + 1)
          IF (I2 .NE. I1) THEN
             D2 = SQRT( REAL(PS + SQ_SUMS(I2) - &
                  2*DOT_PRODUCT(POINT(:),TREE(:,I2)), REAL64) )
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
       ! on which one is closer to the approximation point.
       INNER_CHILD_CLOSER : IF (D1 .LE. D2) THEN
          ! Search the inner child if it could contain a nearer point.
          SEARCH_INNER1 : IF ((F .LT. K) .OR. (D .GT. D1 - RADII(I1))) THEN
             CALL PT_NEAREST_I8(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(2:MID), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_INNER1
          ! Search the outer child if it could contain a nearer point.
          SEARCH_OUTER1 : IF (((F .LT. K) .OR. (D .GT. D2 - RADII(I2))) &
               .AND. (I2 .NE. I1) .AND. (I2 .GT. 0)) THEN
             CALL PT_NEAREST_I8(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(MID+1:), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_OUTER1
       ELSE
          ! Search the outer child if it could contain a nearer point.
          SEARCH_OUTER2 : IF (((F .LT. K) .OR. (D .GT. D2 - RADII(I2))) &
               .AND. (I2 .NE. I1) .AND. (I2 .GT. 0)) THEN
             CALL PT_NEAREST_I8(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(MID+1:), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_OUTER2
          ! Search the inner child if it could contain a nearer point.
          SEARCH_INNER2 : IF ((F .LT. K) .OR. (D .GT. D1 - RADII(I1))) THEN
             CALL PT_NEAREST_I8(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(2:MID), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_INNER2
       END IF INNER_CHILD_CLOSER

    ! Since this is a leaf node, we measure distance to all children.
    ELSE
       DIST_TO_CHILDREN : DO I = 2, SIZE(ORDER)
          ! Measure distance to all children of this node.
          D = SQRT( REAL(PS + SQ_SUMS(ORDER(I)) - &
               2*DOT_PRODUCT(POINT(:),TREE(:,ORDER(I))), REAL64) )
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
       ! This is not the root, We need to pass the updated value of
       ! FOUND back up the recrusion stack.
       FOUND = F
    ELSE
       ! This is the root, initial caller. Sort the distances for return.
       CALL ARGSELECT_R64(DISTS(:F), INDICES(:F), K)
       CALL ARGSORT(DISTS(:K), INDICES(:K))
    END IF SORT_K
  END SUBROUTINE PT_NEAREST_I8


  ! Re-organize a built tree so that it is more usefully packed in memory.
  SUBROUTINE FIX_ORDER_R64(POINTS, SQ_SUMS, RADII, ORDER)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:,:) :: POINTS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:)   :: ORDER
    INTEGER(KIND=INT64) :: I
    ! Execute the swap operation.
    POINTS(:,:) = POINTS(:,ORDER)
    SQ_SUMS(:) = SQ_SUMS(ORDER)
    RADII(:) = RADII(ORDER)
    FORALL (I=1:SIZE(ORDER)) ORDER(I) = I
  END SUBROUTINE FIX_ORDER_R64

  ! Re-organize a built tree so that it is more usefully packed in memory.
  SUBROUTINE FIX_ORDER_I8(POINTS, SQ_SUMS, RADII, ORDER)
    INTEGER(KIND=INT8),  INTENT(INOUT), DIMENSION(:,:) :: POINTS
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:)   :: ORDER
    INTEGER(KIND=INT64) :: I
    ! Execute the swap operation.
    POINTS(:,:) = POINTS(:,ORDER)
    SQ_SUMS(:) = SQ_SUMS(ORDER)
    RADII(:) = RADII(ORDER)
    FORALL (I=1:SIZE(ORDER)) ORDER(I) = I
  END SUBROUTINE FIX_ORDER_I8


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

END MODULE BALL_TREE

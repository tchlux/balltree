MODULE BALL_TREE_R64
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  USE PRUNE,       ONLY: LEVEL
  USE SWAP,        ONLY: SWAP_R64, SWAP_I64
  USE FAST_SELECT, ONLY: ARGSELECT_R64
  USE FAST_SORT,   ONLY: ARGSORT
  IMPLICIT NONE

  ! Max bytes for which a doubling of memory footprint (during copy)
  !  is allowed to happen (switches to using scratch file instead).
  INTEGER(KIND=INT64) :: MAX_COPY_BYTES = 2_INT64 ** 31_INT64

CONTAINS

  ! Re-arrange elements of POINTS into a binary ball tree.
  RECURSIVE SUBROUTINE BUILD_TREE(POINTS, SQ_SUMS, RADII, ORDER,&
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
    REAL(KIND=REAL64) :: MAX_SQ_DIST, SQ_DIST, SHIFT
    EXTERNAL :: DGEMM
    ! Set the leaf size to 1 by default (most possible work required,
    ! but guarantees successful use with any leaf size).
    IF (PRESENT(LEAF_SIZE)) THEN ; LS = LEAF_SIZE
    ELSE                         ; LS = 1
    END IF
    ! If no squared sums were provided, compute them.
    IF (.NOT. PRESENT(COMPUTED_SQ_SUMS) .OR. &
         .NOT. COMPUTED_SQ_SUMS) THEN
       !$OMP PARALLEL DO
       DO I = 1, SIZE(POINTS,2)
          SQ_SUMS(I) = SUM(POINTS(:,I)**2)
       END DO
       !$OMP END PARALLEL DO
    END IF
    ! Set the index of the 'root' of the tree.
    IF (PRESENT(ROOT)) THEN ; CENTER_IDX = ROOT
    ELSE
       ! 1) Compute distances between first point (random) and all others.
       ! 2) Pick the furthest point (on conv hull) from first as the center node.
       J = ORDER(1)
       PT(:) = POINTS(:,J)
       SQ_DISTS(1) = 0.0_REAL64
       !$OMP PARALLEL DO
       ROOT_TO_ALL : DO I = 2, SIZE(ORDER)
          SQ_DISTS(I) = SQ_SUMS(J) + SQ_SUMS(ORDER(I)) - &
               2 * DOT_PRODUCT(POINTS(:,ORDER(I)), PT(:))
       END DO ROOT_TO_ALL
       !$OMP END PARALLEL DO
       CENTER_IDX = MAXLOC(SQ_DISTS(:),1)
       ! Now CENTER_IDX is the selected center for this node in tree.
    END IF

    ! Move the "center" to the first position.
    CALL SWAP_I64(ORDER(1), ORDER(CENTER_IDX))
    ! Measure squared distance beween "center" node and all other points.
    J = ORDER(1)
    PT(:) = POINTS(:,J)
    SQ_DISTS(1) = 0.0_REAL64

    !$OMP PARALLEL DO
    CENTER_TO_ALL : DO I = 2, SIZE(ORDER)
       SQ_DISTS(I) = SQ_SUMS(J) + SQ_SUMS(ORDER(I)) - &
            2 * DOT_PRODUCT(POINTS(:,ORDER(I)), PT(:))
    END DO CENTER_TO_ALL
    !$OMP END PARALLEL DO

    ! Base case for recursion, once we have few enough points, exit.
    IF (SIZE(ORDER) .LE. LS) THEN
       RADII(ORDER(1)) = SQRT(MAXVAL(SQ_DISTS,1))
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
    ! Move the median point (furthest "interior") to the front (inner root).
    CALL SWAP_I64(ORDER(2), ORDER(MID))
    ! Move the furthest point into the spot after the median (outer root).
    CALL SWAP_I64(ORDER(MID+1), ORDER(I))

    !$OMP PARALLEL NUM_THREADS(2)
    !$OMP SECTIONS
    !$OMP SECTION
    ! Recurisively create this tree.
    !   build a tree with the root being the furthest from this center
    !   for the remaining "interior" points of this center node.
    CALL BUILD_TREE(POINTS, SQ_SUMS, RADII, ORDER(2:MID), 1_INT64, LS, .TRUE.)
    !$OMP SECTION
    !   build a tree with the root being the furthest from this center
    !   for the remaining "exterior" points of this center node.
    !   Only perform this operation if there are >0 points available.
    IF (MID < SIZE(ORDER)) &
         CALL BUILD_TREE(POINTS, SQ_SUMS, RADII, &
         ORDER(MID+1:), 1_INT64, LS, .TRUE.)
    !$OMP END SECTIONS
    !$OMP END PARALLEL
  END SUBROUTINE BUILD_TREE


  ! Compute the K nearest elements of TREE to each point in POINTS.
  SUBROUTINE NEAREST(POINTS, K, TREE, SQ_SUMS, RADII, ORDER, &
       LEAF_SIZE, INDICES, DISTS, TO_SEARCH)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: POINTS, TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: DISTS
    INTEGER(KIND=INT64), INTENT(IN),  OPTIONAL    :: TO_SEARCH
    ! Local variables.
    INTEGER(KIND=INT64) :: I, B, BUDGET
    IF (PRESENT(TO_SEARCH)) THEN ; BUDGET = MAX(K, TO_SEARCH)
    ELSE                         ; BUDGET = SIZE(ORDER)       ; END IF
    ! For each point in this set, use the recursive branching
    ! algorithm to identify the nearest elements of TREE.
    !$OMP PARALLEL DO PRIVATE(B)
    DO I = 1, SIZE(POINTS,2)
       B = BUDGET
       CALL PT_NEAREST(POINTS(:,I), K, TREE, SQ_SUMS, RADII, ORDER, &
            LEAF_SIZE, INDICES(:,I), DISTS(:,I), B)
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE NEAREST

  ! Compute the K nearest elements of TREE to each point in POINTS.
  SUBROUTINE PT_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, &
       ORDER, LEAF_SIZE, INDICES, DISTS, CHECKS)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: POINT
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(:) :: DISTS
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: CHECKS
    ! Local variables
    INTEGER(KIND=INT64), DIMENSION(SIZE(ORDER)/2) :: MEASURED_INDS
    INTEGER(KIND=INT64), DIMENSION(SIZE(ORDER)/2) :: MEASURED_SIZES
    REAL(KIND=REAL64), DIMENSION(SIZE(ORDER)/2) :: NEGATED_LOWERS
    REAL(KIND=REAL64), DIMENSION(SIZE(ORDER)/2) :: MEASURED_LOWERS
    INTEGER(KIND=INT64) :: I, S, I1, I2, MID
    INTEGER(KIND=INT64) :: ALLOWED_CHECKS, NEXT_BEST, NUM_MEASURED
    REAL(KIND=REAL64)   :: D, D1, D2, PS
    ! Set the default value for the number of allowed distance calculations.
    IF (PRESENT(CHECKS)) THEN
       ALLOWED_CHECKS = CHECKS
    ELSE
       ALLOWED_CHECKS = SIZE(ORDER)
    END IF
    ! Initialize the indices and distances to all be "null" values.
    INDICES(:) = 0
    DISTS(:) = HUGE(DISTS(1))
    ! Initialize the MEASURED_* arrays for tracking a heap of
    !  subtrees that need to be searched for the nearest neighbors.
    NUM_MEASURED = 0
    MEASURED_INDS(:) = 0
    MEASURED_SIZES(:) = 0
    MEASURED_LOWERS(:) = HUGE(MEASURED_LOWERS(1))
    ! Compute the square sum of the search point
    !  (to accelerate later distance calculations).
    PS = SUM(POINT(:)**2)
    ! Measure the distance to the root node and store it.
    I = ORDER(1)
    D = SQRT(PS + SQ_SUMS(I) - &
         2*DOT_PRODUCT(POINT(:), TREE(:,I)))
    CALL INSERT_MAX(I, D)
    ALLOWED_CHECKS = ALLOWED_CHECKS - 1
    ! Store the statistics about the root node.
    I = 1; S = SIZE(ORDER); CALL INSERT_MIN(I, S, D)
    ! Loop until there are no more potential avenues to descend.
    SEARCH_TREE : DO WHILE (NUM_MEASURED .GT. 0)
       ! Get the next best subtree to check.
       CALL POP_MIN(I, S, D)
       ! If the lower bound for the nearest point in this subtree is
       !  greater than the furthest point found so far, skip this tree.
       IF (D .GE. DISTS(1)) CYCLE SEARCH_TREE
       ! Measure distance to all children if the next best is a leaf.
       IF (S .LE. LEAF_SIZE) THEN
          I1 = I + 1
          I2 = I + S - 1
          DIST_TO_CHILDREN : DO I = I1, I2
             ! Measure distance to all children of this node.
             D = SQRT(PS + SQ_SUMS(ORDER(I)) - &
                  2*DOT_PRODUCT(POINT(:),TREE(:,ORDER(I))))
             ! Store this distance.
             CALL INSERT_MAX(ORDER(I), D)
             ALLOWED_CHECKS = ALLOWED_CHECKS - 1
             IF (ALLOWED_CHECKS .LE. 0) EXIT SEARCH_TREE
          END DO DIST_TO_CHILDREN
       ! Measure distance to the two children and compute bounds.
       ELSE
          ! Measure the distance to the two children and compute the
          !  bounds for nearest possible distance to their descendants.
          I = I + 1 ! <- Get inner child index given root index.
          ! Compute the index of the last point that belongs to the
          !  inner child (and the size of the current set of points).
          MID = (S + 2) / 2
          ! Measure distance to inner child and store.
          I1 = ORDER(I)
          D1 = SQRT(PS + SQ_SUMS(I1) - &
               2*DOT_PRODUCT(POINT(:),TREE(:,I1)))
          CALL INSERT_MAX(I1, D1)
          ALLOWED_CHECKS = ALLOWED_CHECKS - 1
          IF (ALLOWED_CHECKS .LE. 0) EXIT SEARCH_TREE
          ! For subtrees that "contain" the search point, encode the
          !  distance to the root in the distance lower bound.
          IF (D1 .LT. RADII(I1)) THEN
             D1 = -1.0_REAL64 / (1.0_REAL64 + D1)
          ELSE
             D1 = D1 - RADII(I1)
          END IF
          ! Store the descriptive variables for the inner child lower bound.
          IF (D1 .LT. DISTS(1)) CALL INSERT_MIN(I, MID-1, D1)
          ! Check to see if there actually is an outer child.
          IF (MID+1 .LE. S) THEN
             I = I + MID - 1 ! <- the index immediately after MID
             ! Measure distance to outer child and store.
             I2 = ORDER(I)
             D2 = SQRT(PS + SQ_SUMS(I2) - &
                  2*DOT_PRODUCT(POINT(:),TREE(:,I2)))
             CALL INSERT_MAX(I2, D2)
             ALLOWED_CHECKS = ALLOWED_CHECKS - 1
             IF (ALLOWED_CHECKS .LE. 0) EXIT SEARCH_TREE
             ! For subtrees that "contain" the search point, encode the
             !  distance to the root in the measurement.
             IF (D2 .LT. RADII(I2)) THEN
                D2 = -1.0_REAL64 / (1.0_REAL64 + D2)
             ELSE
                D2 = D2 - RADII(I2)
             END IF
             ! Store the descriptive variables for the outer child.
             IF (D2 .LT. DISTS(1)) CALL INSERT_MIN(I, S-MID, D2)
          END IF
       END IF
    END DO SEARCH_TREE
    ! Sort the distances for return (they were in a heap before).
    CALL ARGSORT(DISTS(:), INDICES(:))

  CONTAINS

    ! Insert a new element into a min-heap.
    SUBROUTINE INSERT_MIN(INDEX, SUBTREE_SIZE, LOWER)
      INTEGER(KIND=INT64), INTENT(IN) :: INDEX, SUBTREE_SIZE
      REAL(KIND=REAL64),   INTENT(IN) :: LOWER
      ! Temparary index.
      INTEGER(KIND=INT64) :: I, IP
      ! Insert at the end.
      NUM_MEASURED = NUM_MEASURED + 1
      MEASURED_INDS(NUM_MEASURED) = INDEX
      MEASURED_SIZES(NUM_MEASURED) = SUBTREE_SIZE
      MEASURED_LOWERS(NUM_MEASURED) = LOWER
      ! Compare upwards to maintain HEAP property.
      I = NUM_MEASURED
      CHECK_MIN_HEAP : DO
         IP = I / 2
         IF (IP .LT. 1) THEN
            EXIT CHECK_MIN_HEAP            
         ELSE IF (MEASURED_LOWERS(IP) .LE. LOWER) THEN
            EXIT CHECK_MIN_HEAP
         ELSE
            CALL SWAP_I64(MEASURED_INDS(IP), MEASURED_INDS(I))
            CALL SWAP_I64(MEASURED_SIZES(IP), MEASURED_SIZES(I))
            CALL SWAP_R64(MEASURED_LOWERS(IP), MEASURED_LOWERS(I))
            I = IP
         END IF
      END DO CHECK_MIN_HEAP
    END SUBROUTINE INSERT_MIN

    ! Pop the root from a min-heap.
    SUBROUTINE POP_MIN(INDEX, SUBTREE_SIZE, LOWER)
      INTEGER(KIND=INT64), INTENT(OUT) :: INDEX, SUBTREE_SIZE
      REAL(KIND=REAL64),   INTENT(OUT) :: LOWER
      ! Temporary index.
      INTEGER(KIND=INT64) :: I, I1, I2
      ! Get the minimum element.
      INDEX = MEASURED_INDS(1)
      SUBTREE_SIZE = MEASURED_SIZES(1)
      LOWER = MEASURED_LOWERS(1)
      ! Swap the last element to the root.
      CALL SWAP_I64(MEASURED_INDS(1), MEASURED_INDS(NUM_MEASURED))
      CALL SWAP_I64(MEASURED_SIZES(1), MEASURED_SIZES(NUM_MEASURED))
      CALL SWAP_R64(MEASURED_LOWERS(1), MEASURED_LOWERS(NUM_MEASURED))
      MEASURED_INDS(NUM_MEASURED) = 0
      MEASURED_SIZES(NUM_MEASURED) = -1
      MEASURED_LOWERS(NUM_MEASURED) = 0.0
      NUM_MEASURED = NUM_MEASURED - 1
      ! Comparex downwards to maintain HEAP property.
      I = 1
      CHECK_MIN_HEAP : DO
         ! Get the index ofthe smaller child into I1.
         I1 = 2*I
         I2 = I1 + 1
         IF (I1 .GT. NUM_MEASURED) THEN
            EXIT CHECK_MIN_HEAP
         ELSE IF (I2 .GT. NUM_MEASURED) THEN
            CONTINUE
         ELSE IF (MEASURED_LOWERS(I2) .LT. MEASURED_LOWERS(I1)) THEN
            I1 = I2
         END IF
         ! Swap down if the smaller child is less than the current.
         IF (MEASURED_LOWERS(I1) .LT. MEASURED_LOWERS(I)) THEN
            CALL SWAP_I64(MEASURED_INDS(I1), MEASURED_INDS(I))
            CALL SWAP_I64(MEASURED_SIZES(I1), MEASURED_SIZES(I))
            CALL SWAP_R64(MEASURED_LOWERS(I1), MEASURED_LOWERS(I))
            I = I1
         ELSE
            EXIT CHECK_MIN_HEAP
         END IF
      END DO CHECK_MIN_HEAP
    END SUBROUTINE POP_MIN

    ! Insert a new element into a min-heap.
    SUBROUTINE INSERT_MAX(INDEX, DISTANCE)
      INTEGER(KIND=INT64), INTENT(IN) :: INDEX
      REAL(KIND=REAL64),   INTENT(IN) :: DISTANCE
      INTEGER(KIND=INT64) :: I, IP, I1, I2
      ! Do not add the new point if its distance is greater.
      IF (DISTANCE .GE. DISTS(1)) RETURN
      ! Otherwise, replace the existing max with the new value.
      INDICES(1) = INDEX
      DISTS(1) = DISTANCE
      ! Comparex downwards to maintain HEAP property.
      I = 1
      CHECK_MAX_HEAP_AFTER_REMOVE : DO
         ! Get the index ofthe smaller child into I1.
         I1 = 2*I
         I2 = I1 + 1
         IF (I1 .GT. K) THEN
            EXIT CHECK_MAX_HEAP_AFTER_REMOVE
         ELSE IF (I2 .GT. K) THEN
            I1 = I2
         ELSE IF (DISTS(I2) .GT. DISTS(I1)) THEN
            I1 = I2
         END IF
         ! Swap down if the larger child is greater than the current.
         IF (DISTS(I1) .GT. DISTS(I)) THEN
            CALL SWAP_I64(INDICES(I1), INDICES(I))
            CALL SWAP_R64(DISTS(I1), DISTS(I))
            I = I1
         ELSE
            EXIT CHECK_MAX_HEAP_AFTER_REMOVE
         END IF
      END DO CHECK_MAX_HEAP_AFTER_REMOVE
    END SUBROUTINE INSERT_MAX

  END SUBROUTINE PT_NEAREST


  ! Re-organize a built tree so that it is more usefully packed in memory.
  SUBROUTINE FIX_ORDER(POINTS, SQ_SUMS, RADII, ORDER, COPY)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:,:) :: POINTS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:)   :: ORDER
    LOGICAL, INTENT(IN), OPTIONAL :: COPY
    LOGICAL :: SHOULD_COPY
    INTEGER(KIND=INT64) :: I
    ! Default to copy (in memory) if there is less than 1 GB of data.
    IF (PRESENT(COPY)) THEN ; SHOULD_COPY = COPY
    ELSE ; SHOULD_COPY = SIZEOF(POINTS) .LE. MAX_COPY_BYTES
    END IF
    ! Reorder all of the data. Use a scratch file for large data sets.
    IF (SHOULD_COPY) THEN
       POINTS(:,:) = POINTS(:,ORDER)
    ELSE
       ! Open scratch file for writing all of the points in order.
       OPEN(UNIT=1, STATUS='SCRATCH', ACTION='READWRITE', FORM='UNFORMATTED', ACCESS='STREAM')
       ! Write all points to a scratch file in the correct order.
       DO I = 1, SIZE(ORDER)
          WRITE(UNIT=1) POINTS(:,ORDER(I))
       END DO
       ! Read all points from file (they are now ordered correctly).
       READ(UNIT=1,POS=1) POINTS(:,:)
       ! Close scratch file.
       CLOSE(UNIT=1)
    END IF
    ! Always copy the square sums and the radii in memory (shouldn't be bad).
    SQ_SUMS(:) = SQ_SUMS(ORDER)
    RADII(:) = RADII(ORDER)
    ! Reset the order because now it is the expected format.
    FORALL (I=1:SIZE(ORDER)) ORDER(I) = I
  END SUBROUTINE FIX_ORDER

  ! ------------------------------------------------------------------
  !               APPROXIMATE NEAREST NEIGHBOR LOOKUP
  ! ------------------------------------------------------------------

  ! Compute the K nearest elements of TREE to each point in POINTS.
  SUBROUTINE APPROX_NEAREST(POINTS, K, TREE, SQ_SUMS, RADII, ORDER, &
       LEAF_SIZE, INDICES, DISTS, LOOK_AHEAD, RANDOMIZED)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: POINTS, TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: DISTS
    INTEGER(KIND=INT64), INTENT(IN),  OPTIONAL    :: LOOK_AHEAD
    LOGICAL, INTENT(IN), OPTIONAL :: RANDOMIZED
    ! Local variables.
    LOGICAL :: RANDOM_TRAVERSAL
    INTEGER(KIND=INT64) :: I, B, LH
    INTEGER(KIND=INT64), DIMENSION(:), ALLOCATABLE :: INDS_BUFFER
    REAL(KIND=REAL64),   DIMENSION(:), ALLOCATABLE :: DISTS_BUFFER
    IF (PRESENT(LOOK_AHEAD)) THEN ; LH = LOOK_AHEAD
    ELSE ; LH = MIN(5, CEILING(LOG(REAL(SIZE(ORDER))) / LOG(2.0))) ; END IF
    IF (PRESENT(RANDOMIZED)) THEN ; RANDOM_TRAVERSAL = RANDOMIZED
    ELSE ; RANDOM_TRAVERSAL = .FALSE. ; END IF
    ! Allocate the buffers for holding the nearest indices and distances.
    ALLOCATE( INDS_BUFFER(1:K+2**LH), DISTS_BUFFER(1:K+2**LH) )
    ! For each point in this set, use the recursive branching
    ! algorithm to identify the nearest elements of TREE.
    !$OMP PARALLEL DO PRIVATE(INDS_BUFFER, DISTS_BUFFER, B)
    DO I = 1, SIZE(POINTS,2)
       CALL PT_APPROX_NEAREST(POINTS(:,I), K, TREE, SQ_SUMS, RADII, &
            ORDER, LEAF_SIZE, INDS_BUFFER, DISTS_BUFFER, LH, &
            RANDOMIZED=RANDOM_TRAVERSAL)
       ! Sort the first K elements of the temporary arry for return.
       INDICES(:,I) = INDS_BUFFER(:K)
       DISTS(:,I) = DISTS_BUFFER(:K)
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE APPROX_NEAREST

  ! Compute the K nearest elements of TREE to each point in POINTS
  ! using the "look ahead" method for determining tree traversal.
  RECURSIVE SUBROUTINE PT_APPROX_NEAREST(POINT, K, TREE, SQ_SUMS, &
       RADII, ORDER, LEAF_SIZE, INDICES, DISTS, LOOK_AHEAD, FOUND, &
       PT_SS, RANDOMIZED)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: POINT
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE, LOOK_AHEAD
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(:) :: DISTS
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: FOUND
    REAL(KIND=REAL64),   INTENT(IN),    OPTIONAL   :: PT_SS
    LOGICAL, INTENT(IN), OPTIONAL :: RANDOMIZED
    ! Local variables
    INTEGER(KIND=INT64) :: F, I, I1, I2, MID, LMID
    REAL(KIND=REAL64)   :: D, D1, D2
    REAL(KIND=REAL64)   :: PS
    LOGICAL :: RANDOM_TRAVERSAL
    ! Storage space for the inner and outer children indices / distances.
    INTEGER(KIND=INT64), DIMENSION(1:2**LOOK_AHEAD) :: LEVEL_INDS
    REAL(KIND=REAL64), DIMENSION(1:2**LOOK_AHEAD) :: LEVEL_DISTS
    ! Initialize randomization.
    IF (PRESENT(RANDOMIZED)) THEN
       RANDOM_TRAVERSAL = RANDOMIZED
    ELSE
       RANDOM_TRAVERSAL = .FALSE.
    END IF
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
    ! If there are enough layers beneath this one, then look ahead.
    KEEP_LOOKING : IF (SIZE(ORDER) .GE. 2**LOOK_AHEAD) THEN
       ! If this is NOT a leaf node, then recurse.
       BRANCH_OR_LEAF : IF (SIZE(ORDER) .GT. 2*LEAF_SIZE+1) THEN
          ! Compute the indices of the inner and outer children.
          MID = (SIZE(ORDER) + 2) / 2
          LMID = 2 ** (LOOK_AHEAD - 1)
          I1 = 0_INT64
          CALL LEVEL(MID-1, LOOK_AHEAD-1, LEVEL_INDS(:LMID), I1, 2_INT64)
          I2 = 0_INT64
          CALL LEVEL(SIZE(ORDER)-MID, LOOK_AHEAD-1, LEVEL_INDS(LMID+1:), I2, MID+1)
          ! Compute distance to the LOOK_AHEAD next level.
          ! Traverse down the branch that looks better.
          ! 
          ! Measure distance to all inner children.
          DO I = 1, I1
             LEVEL_DISTS(I) = SQRT(PS + SQ_SUMS(ORDER(LEVEL_INDS(I))) - &
                  2*DOT_PRODUCT(POINT(:),TREE(:,ORDER(LEVEL_INDS(I)))))
             ! Store this measured distance.
             F = F + 1
             INDICES(F) = ORDER(LEVEL_INDS(I))
             DISTS(F) = LEVEL_DISTS(I)
          END DO
          ! Measure distance to all outer children.
          DO I = LMID+1, LMID+I2
             LEVEL_DISTS(I) = SQRT(PS + SQ_SUMS(ORDER(LEVEL_INDS(I))) - &
                  2*DOT_PRODUCT(POINT(:),TREE(:,ORDER(LEVEL_INDS(I)))))
             ! Store this measured distance.
             F = F + 1
             INDICES(F) = ORDER(LEVEL_INDS(I))
             DISTS(F) = LEVEL_DISTS(I)
          END DO
          ! Reduce the kept points to only those that are closest.
          CALL ARGSELECT_R64(DISTS(:F), INDICES(:F), MIN(K,F))
          F = MIN(K, F)
          ! Compute the probability that the inner (left) is less than the outer (right).
          CALL PROBABILITY_LESS_THAN(LEVEL_DISTS(1:I1), LEVEL_DISTS(LMID+1:LMID+I2), D1)
          ! Compute the probability that the outer (right) is less than the inner (left).
          CALL PROBABILITY_LESS_THAN(LEVEL_DISTS(LMID+1:LMID+I2), LEVEL_DISTS(1:I1), D2)
          ! Compute the probability of traversing for randomized search.
          IF (RANDOM_TRAVERSAL) THEN
             ! Normalize the probability that D1 is less.
             D1 = D1 / (D1 + D2)
             ! Generate a random number to decide which path to traverse.
             CALL RANDOM_NUMBER(D)
             ! Determine whether or not D1 is less than D2 based on randomness.
             IF (D1 .LE. D) THEN
                D1 = 2.0_REAL64
             ELSE
                D1 = -1.0_REAL64
             END IF
          END IF
          ! Traverse inner (left).
          IF (D1 .GE. D2) THEN
             CALL PT_APPROX_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, ORDER(2:MID), &
                  LEAF_SIZE, INDICES, DISTS, LOOK_AHEAD, F, PS, RANDOM_TRAVERSAL)
          ! Traverse outer (right).
          ELSE
             CALL PT_APPROX_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, ORDER(MID+1:), &
                  LEAF_SIZE, INDICES, DISTS, LOOK_AHEAD, F, PS, RANDOM_TRAVERSAL)
          END IF
       ELSE
          ! This is a leaf, distances need to be measured to all the
          ! children of each child node.
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
    ELSE
       ! There are no further levels to measure distance to, exit.
       RETURN
    END IF KEEP_LOOKING

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
  END SUBROUTINE PT_APPROX_NEAREST

  ! Given two lists of numbers, compute the probability that a random
  !  number selected from the first list is less than a random number
  !  selected from the second list.
  SUBROUTINE PROBABILITY_LESS_THAN(L1, L2, P)
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:) :: L1, L2
    REAL(KIND=REAL64), INTENT(OUT) :: P
    ! Local variables.
    INTEGER(KIND=INT64), DIMENSION(SIZE(L2)) :: INDS
    INTEGER(KIND=INT64) :: I, J, COUNT
    ! Sort the second list.
    CALL ARGSORT(L2(:), INDS(:))
    ! Sum the probabilities that each element in L1 is less than a
    !  a random element for L2, then divide by the size of L1.
    COUNT = 0
    DO I = 1, SIZE(L1)
       CALL FIRST_GREATER_OR_EQUAL(L2, L1(I), J)
       COUNT = COUNT + (SIZE(L2) - J + 1)
    END DO
    P = REAL(COUNT, KIND=REAL64) / REAL(SIZE(L2)*SIZE(L1), KIND=REAL64)

  CONTAINS
    ! Identify the index of the first element in a list that is greater
    !  than or equal to a given value.
    SUBROUTINE FIRST_GREATER_OR_EQUAL(LIST, VAL, INDEX)
      REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: LIST
      REAL(KIND=REAL64), INTENT(IN) :: VAL
      INTEGER(KIND=INT64), INTENT(OUT) :: INDEX
      ! Local variables.
      INTEGER(KIND=INT64) :: FRONT, BACK
      ! Check for special case where every value in the list is greater.
      IF (LIST(1) .GE. VAL) THEN
         INDEX = 1
         RETURN
      END IF
      ! Check for special case where every value in the list is less.
      IF (LIST(SIZE(LIST)) .LT. VAL) THEN
         INDEX = SIZE(LIST)+1
         RETURN
      END if
      ! Perform a binary search to find the first index of an element
      ! in the list that is strictly less than the given value.
      FRONT = 1
      BACK = SIZE(L2)
      binary_search : DO WHILE (FRONT .LT. BACK)
         INDEX = (FRONT + BACK) / 2
         ! Check loop termination condition.
         IF ((INDEX .EQ. FRONT) .OR. (INDEX .EQ. BACK)) THEN
            IF (LIST(FRONT) .GE. VAL) THEN
               INDEX = FRONT
            ELSE
               INDEX = BACK
            END IF
            RETURN
         ! Shrink front and back towards each other.
         ELSE
            IF (LIST(INDEX) .LT. VAL) THEN
               FRONT = INDEX
            ELSE
               BACK = INDEX
            END IF
         END IF
      END DO binary_search
    END SUBROUTINE FIRST_GREATER_OR_EQUAL
  END SUBROUTINE PROBABILITY_LESS_THAN


END MODULE BALL_TREE_R64

# TODO: Add code for incrementally adding points to the tree
# TODO: Add indexing scheme (after building the tree? links through layers?)
# TODO: Add a search that takes advantage of indexing.

# --------------------------------------------------------------------
#          DECLARING DIFFERENT TYPE BACK ENDS FOR BALL TREE
# 
# # Declare the methods based on the dtype.
# if ('int8' in str(self.ttype)):
#     if ('uint8' in str(self.ttype)):
#         import warnings
#         warnings.warn("This ball tree only handles signed integers. Make sure to subtract 128 from all provided values before using this code.")
#     self.sstype = np.int64
#     self._build_tree = ball_tree_i8.build_tree
#     self._fix_order  = ball_tree_i8.fix_order
#     self._bt_nearest = ball_tree_i8.nearest
#     # TODO: Need to write "approx nearest" function for I8.
#     self._bt_approx_nearest = lambda *args, **kwargs: print("ERROR: Unsupported operation.")
# --------------------------------------------------------------------


# --------------------------------------------------------------------------
# Alternative partitioning strategies tested, not better enough to merit use
#  N = 100k
#  D = 100
#  K = 5
# 
# # Move a random interior point into the root location for the inner set.
# CALL RANDOM_NUMBER(SHIFT)
# I = SHIFT * (MID-2) + 2
# CALL SWAP_I64(ORDER(2), ORDER(I))
# # Move a random exterior point into the root location for the outer set.
# CALL RANDOM_NUMBER(SHIFT)
# I = MID + (SHIFT * (SIZE(ORDER)-MID-1) + 1)
# CALL SWAP_I64(ORDER(MID+1), ORDER(I))
# # .022 .020 .021 .020 .029 .020

# # Move the closest interior point into the root location for the inner set.
# I = MINLOC(SQ_DISTS(2:MID),1) + 1
# CALL SWAP_I64(ORDER(2), ORDER(I))
# # Move the closest exterior point into the root location for the outer set.
# I = MID + MINLOC(SQ_DISTS(MID+1:),1)
# CALL SWAP_I64(ORDER(MID+1), ORDER(I))
# # .02 .02 .02 .02 .02
# --------------------------------------------------------------------------


# -----------------------------------------------------------------------------
#      OLD FORTRAN METHOD FOR DOING APPROXIMATE NEAREST NEIGHBOR
# 
# ! Find the minimum of the two.
# IF (MINVAL(LEVEL_DISTS(1:I1)) .LE. MINVAL(LEVEL_DISTS(LMID+1:LMID+I2))) THEN
#    ! The left is better.
#    CALL PT_APPROX_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, ORDER(2:MID), &
#         LEAF_SIZE, INDICES, DISTS, LOOK_AHEAD, FOUND, PT_SS, RANDOM_TRAVERSAL)
# ELSE
#    ! The right is better.
#    CALL PT_APPROX_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, ORDER(MID+1:), &
#         LEAF_SIZE, INDICES, DISTS, LOOK_AHEAD, FOUND, PT_SS, RANDOM_TRAVERSAL)
# END IF
# -----------------------------------------------------------------------------


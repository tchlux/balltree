# TODO: Add code for incrementally adding points to the tree
# TODO: Add indexing scheme (after building the tree? links through layers?)
# TODO: Add a search that takes advantage of indexing.

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

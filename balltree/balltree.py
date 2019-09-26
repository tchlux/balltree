# Ball Tree wrapper for Fortran module.

# TODO: Add another operation mode for a 'byte' tree.
# TODO: Mix test cases for different type trees, make one testing suite.
# TODO: If data is 'uint8' then add 128 to the values before building tree.
# TODO: Add code for incrementally adding points to the tree

import os 
import fmodpy
import numpy as np
PATH_TO_HERE   = os.path.dirname(os.path.abspath(__file__))
# Import the Fortran utilities.
PATH_TO_BT     = os.path.join(PATH_TO_HERE, "ball_tree.f90")
ball_tree = fmodpy.fimport(PATH_TO_BT, output_directory=PATH_TO_HERE,
                           autocompile_extra_files=True)
PATH_TO_SORT   = os.path.join(PATH_TO_HERE, "fast_sort.f90")
fast_sort = fmodpy.fimport(PATH_TO_SORT, output_directory=PATH_TO_HERE,
                           autocompile_extra_files=True)
PATH_TO_SELECT = os.path.join(PATH_TO_HERE, "fast_select.f90")
fast_select = fmodpy.fimport(PATH_TO_SELECT, output_directory=PATH_TO_HERE,
                             autocompile_extra_files=True)


# ------------------------------------------------------------------
#                        FastSort method
# 
# This routine uses a combination of QuickSort (with modestly
# intelligent pivot selection) and Insertion Sort (for small arrays)
# to achieve very fast average case sort times for both random and
# partially sorted data. The pivot is selected for QuickSort as the
# median of the first, middle, and last values in the array.
# 
# Arguments:
# 
#   VALUES   --  A 1D array of real numbers.
# 
# Optional:
# 
#   INDICES  --  A 1D array of original indices for elements of VALUES.
#   MIN_SIZE --  An positive integer that represents the largest
#                sized VALUES for which a partition about a pivot
#                is used to reduce the size of a an unsorted array.
#                Any size less than this will result in the use of
#                INSERTION_ARGSORT instead of ARGPARTITION.
# 
# Output:
# 
#   The elements of the array VALUES are sorted and all elements of
#   INDICES are sorted symmetrically (given INDICES = 1, ...,
#   SIZE(VALUES) beforehand, final INDICES will show original index
#   of each element of VALUES before the sort operation).
# 
def argsort(values, indices=None, min_size=None):
    indices = np.arange(len(values))
    return fast_sort.argsort(values, indices, min_size=min_size)


# ------------------------------------------------------------------
#                       FastSelect method
# 
# Given VALUES list of numbers, rearrange the elements of VALUES
# such that the element at index K has rank K (holds its same
# location as if all of VALUES were sorted). Symmetrically rearrange
# array INDICES to keep track of prior indices.
# 
# This algorithm uses the same conceptual approach as Floyd-Rivest,
# but instead of standard-deviation based selection of bounds for
# recursion, a rank-based method is used to pick the subset of
# values that is searched. This simplifies the code and improves
# readability, while achieving the same tunable performance.
# 
# Arguments:
# 
#   VALUES   --  A 1D array of real numbers.
#   K        --  A positive integer for the rank index about which
#                VALUES should be rearranged.
# Optional:
# 
#   INDICES  --  A 1D array of original indices for elements of VALUES.
#   DIVISOR  --  A positive integer >= 2 that represents the
#                division factor used for large VALUES arrays.
#   MAX_SIZE --  An integer >= DIVISOR that represents the largest
#                sized VALUES for which the worst-case pivot value
#                selection is tolerable. A worst-case pivot causes
#                O( SIZE(VALUES)^2 ) runtime. This value should be
#                determined heuristically based on compute hardware.
# 
# Output:
# 
#   The elements of the array VALUES are rearranged such that the
#   element at position VALUES(K) is in the same location it would
#   be if all of VALUES were in sorted order. Also known as,
#   VALUES(K) has rank K.
# 
def argselect(values, k, indices=None, divisor=None, max_size=None):
    if (indices is None): indices = np.arange(len(values))
    return fast_select.argselect_r64(values, indices, k+1, divisor=divisor, max_size=max_size)


# Class for constructing a ball tree.
class BallTree:
    # Given points and a leaf size, construct a ball tree.
    def __init__(self, points, leaf_size=1, transpose=True):
        if transpose:
            points = points.T
        # Handle different use cases.
        if ('int8' in points.dtype.name):
            self.leaf_size = leaf_size
            self.tree    = np.asarray(points, order='F', dtype='int8')
            self.sq_sums = np.zeros(self.tree.shape[1],  dtype='int64')
            self.order   = np.arange(self.tree.shape[1], dtype='int64') + 1
            self.build_tree = ball_tree.build_tree_i8
            self.fix_order  = ball_tree.fix_order_i8
            self.bt_nearest = ball_tree.nearest_i8
            self.dtype = np.int8
            self.itype = np.int64
        else:
            self.leaf_size = leaf_size
            self.tree    = np.asarray(points, order='F', dtype='float64')
            self.sq_sums = np.zeros(self.tree.shape[1],  dtype='float64')
            self.order   = np.arange(self.tree.shape[1], dtype='int64') + 1
            self.build_tree = ball_tree.build_tree_r64
            self.fix_order  = ball_tree.fix_order_r64
            self.bt_nearest = ball_tree.nearest_r64
            self.dtype = np.float64
            self.itype = np.int64
        self.radii   = np.zeros(self.tree.shape[1],  dtype='float64')
        # Build tree (in-place operation).
        #    .tree     will not be modified
        #    .sq_sums  will contain the squared sums of each point
        #    .radii    will be modified to have the radius of specific
        #              node, or 0.0 if this point is a leaf.
        #    .order    will be the list of indices (1-indexed) that
        #              determine the structure of the ball tree.
        self.build_tree(self.tree, self.sq_sums, self.radii,
                        self.order, leaf_size=self.leaf_size)
        self.index_mapping = self.order.copy()-1
        # Restructure the ball tree so the points are in locally
        # contiguous blocks of memory (local by branch + leaf).
        self.fix_order(self.tree, self.sq_sums, self.radii, self.order)

    # Find the "k" nearest neighbors to all points in z. Uses the same
    # interface as the "BallTree.query" function, see help for more info.
    def nearest(self, *args, **kwargs): return self.query(*args, **kwargs)

    # Find the "k" nearest neighbors to all points in z.
    def query(self, z, k=1, leaf_size=None, return_distance=True, transpose=True):
        if (len(z.shape) != 2) or (z.shape[1] != self.tree.shape[0]):
            class WrongDimension(Exception): pass
            raise(WrongDimension(f"The provided query points must be in a row-vector matrix with dimension (..,{self.tree.shape[0]})."))
        if (leaf_size is None): leaf_size = self.leaf_size
        if transpose: z = z.T
        k = min(k, self.tree.shape[1])
        # Initialize holders for output.
        points  = np.asarray(z, order='F', dtype=self.dtype)
        indices = np.ones((k, points.shape[1]), order='F', dtype=self.itype)
        dists   = np.ones((k, points.shape[1]), order='F', dtype='float64')
        # Compute the nearest neighbors.
        self.bt_nearest(points, k, self.tree, self.sq_sums, self.radii,
                        self.order, leaf_size, indices, dists)
        # Return the results.
        if return_distance:
            if transpose: return dists.T, indices.T - 1
            else:         return dists,   indices - 1
        else:
            if transpose: return indices.T - 1
            else:         return indices - 1

# Ball Tree wrapper for Fortran module.

# TODO: Add another operation mode for a 'byte' tree.
# TODO: Mix test cases for different type trees, make one testing suite.
# TODO: If data is 'uint8' then add 128 to the values before building tree.
# TODO: Add code for incrementally adding points to the tree

import os 
import fmodpy
import numpy as np

PATH_TO_HERE = os.path.dirname(os.path.abspath(__file__))

# Allow OpenMP to create nested threads.
from multiprocessing import cpu_count
# Set OpenMP variables to allow for greatest parallelism when building tree.
if "OMP_NUM_THREADS" not in os.environ:
    os.environ["OMP_NUM_THREADS"] = str(cpu_count())
if "OMP_MAX_ACTIVE_LEVELS" not in os.environ:
    os.environ["OMP_MAX_ACTIVE_LEVELS"] = str(int(np.ceil(np.log2(cpu_count()))))
if "OMP_NESTED" not in os.environ:
    os.environ["OMP_NESTED"] = "TRUE"

# Define paths to the Fortran utilities.
R64_DIR = os.path.join(PATH_TO_HERE, "r64i64")
PATH_TO_BT_R64 = os.path.join(R64_DIR, "ball_tree.f90")
PATH_TO_PRUNE = os.path.join(R64_DIR, "prune.f90")
PATH_TO_SORT   = os.path.join(R64_DIR, "fast_sort.f90")
PATH_TO_SELECT = os.path.join(R64_DIR, "fast_select.f90")
# PATH_TO_BT_I8 = os.path.join(PATH_TO_HERE, "ball_tree_i8.f90")

# Import the Fortran utilities.
dependencies = ["swap.f90", "prune.f90", "fast_select.f90", "fast_sort.f90", "ball_tree.f90"]
ball_tree_r64  = fmodpy.fimport(PATH_TO_BT_R64, output_dir=R64_DIR, omp=True,
                                verbose=False, depends_files=dependencies,
                                autocompile=True).ball_tree_r64
prune = fmodpy.fimport(PATH_TO_PRUNE, output_dir=R64_DIR, verbose=False).prune
fast_sort = fmodpy.fimport(PATH_TO_SORT, output_dir=R64_DIR, verbose=False).fast_sort
fast_select = fmodpy.fimport(PATH_TO_SELECT, output_dir=R64_DIR, verbose=False).fast_select
# ball_tree_i8  = fmodpy.fimport(PATH_TO_BT_I8, output_dir=PATH_TO_HERE, omp=True)

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
    def __init__(self, points=None, dtype=None, transpose=True,
                 leaf_size=1, reorder=True, root=None):
        self.leaf_size = leaf_size
        # Assign the data type if it was provided.
        if (dtype is not None):
            self.ttype = dtype
            self._set_type_internals()
        else: self.ttype = None
        # Assign the various pruning functions.
        self._inner = prune.inner
        self._outer = prune.outer
        self._top = prune.top
        # Assign the points and build the tree, if provided.
        if (points is not None):
            self.add(points, transpose=transpose)
            self.build(reorder=reorder, root=root)

    # Given the data type of this class, setup which internal methods
    # are called (for the appropriate number sizes).
    def _set_type_internals(self):
        # Declare the methods based on the dtype.
        if ('int8' in str(self.ttype)):
            if ('uint8' in str(self.ttype)):
                import warnings
                warnings.warn("This ball tree only handles signed integers. Make sure to subtract 128 from all provided values before using this code.")
            self.sstype = np.int64
            self._build_tree = ball_tree_i8.build_tree
            self._fix_order  = ball_tree_i8.fix_order
            self._bt_nearest = ball_tree_i8.nearest
            # TODO: Need to write "approx nearest" function for I8.
            self._bt_approx_nearest = lambda *args, **kwargs: print("ERROR: Unsupported operation.")
        elif ('float64' in str(self.ttype)):
            self.sstype = np.float64
            self._build_tree = ball_tree_r64.build_tree
            self._fix_order  = ball_tree_r64.fix_order
            self._bt_nearest = ball_tree_r64.nearest
            self._bt_approx_nearest = ball_tree_r64.approx_nearest
        else:
            class UnsupportedType(Exception): pass
            raise(UnsupportedType(f"The type '{dtype.name}' is not supported."))

    # Based on the size and 'built' amount and the 'size' internal: if
    # searching, return 'built' portion of tree; if building return
    # the 'size' portion of tree.
    def _get_order(self, search=False):
        # For a search, only use the built portion of the tree.
        if search:
            if (self.built == self.tree.shape[1]): return self.order
            else: return self.order[:self.built]
        # For a build, only use the 'size' kept points in the tree.
        else:
            if (self.size == self.tree.shape[1]): return self.order
            else: return self.order[:self.size]

    # Add points to this ball tree, if the type is not yet defined,
    # then initialize the type of this tree to be same as added points.
    def add(self, points, transpose=True):
        # Assign the data type if it is not set.
        if (self.ttype is None):
            if ('int8' in points.dtype.name): self.ttype = np.int8
            else:                             self.ttype = np.float64
            self._set_type_internals()
        elif (points.dtype != self.ttype):
            class WrongType(Exception): pass
            raise(WrongType(f"The points provided were type '{points.dtype.name}', expected '{self.ttype.name}'."))
        # Transpose the points if the are not already column-vectors.
        if transpose: points = points.T
        # If there are existing points, make sure the new ones match.
        if hasattr(self, "size"):
            # If possible, pack the points into available space.
            if (points.shape[1] <= self.tree.shape[1] - self.size):
                self.tree[:,self.size:self.size+point.shape[1]] = points
                self.size += points.shape[1]
            else:
                old_tree = self.tree
                old_points = self.tree[:,:self.size]
                dim = self.tree.shape[1]
                self.size += points.shape[1]
                self.tree = np.zeros((dim, self.size), dtype=self.ttype)
                # Assign the relevant internals for this tree.
                self.sq_sums = np.zeros(self.tree.shape[1],  dtype=self.sstype)
                self.order   = np.arange(self.tree.shape[1], dtype='int64') + 1
                self.radii   = np.zeros(self.tree.shape[1],  dtype='float64')
                # Pack the old points and the new points into a single tree.
                self.tree[:,:old_points.shape[1]] = old_points
                self.tree[:,old_points.shape[1]:] = points
                # Delete the old tree.
                del old_tree
                self.built = 0
                return
        else:
            # Save the points internally as the tree.
            self.tree    = np.asarray(points, order='F', dtype=self.ttype)
            # Assign the relevant internals for this tree.
            self.sq_sums = np.zeros(self.tree.shape[1],  dtype=self.sstype)
            self.order   = np.arange(self.tree.shape[1], dtype='int64') + 1
            self.radii   = np.zeros(self.tree.shape[1],  dtype='float64')
            # Store BallTree internals for knowing how to evaluate.
            self.built = 0
            self.size = self.tree.shape[1]

    # Build a tree out.
    def build(self, leaf_size=None, reorder=True, has_sq_sums=False, root=None):
        # Get the leaf size if it was not given.
        if (leaf_size is None): leaf_size = self.leaf_size
        # Translate the root from python index to fortran index if provided.
        if (root is not None): root += 1
        # Automatically set 'reorder' to False if the array is more than
        # 1 GB in size because this requires a doubling of memory usage.
        if (reorder is None):   reorder = (self.tree.nbytes / 2**30) <= 1
        # Build tree (in-place operation).
        #    .tree     will not be modified
        #    .sq_sums  will contain the squared sums of each point
        #    .radii    will be modified to have the radius of specific
        #              node, or 0.0 if this point is a leaf.
        #    .order    will be the list of indices (1-indexed) that
        #              determine the structure of the ball tree.
        order = self._get_order(search=False)
        self._build_tree(self.tree, self.sq_sums, self.radii, order,
                         leaf_size=leaf_size, computed_sq_sums=has_sq_sums, 
                         root=root)
        self.built = len(order)
        # Store the index mapping from the build.
        if hasattr(self, "index_mapping"): del self.index_mapping
        self.index_mapping = order.copy() - 1
        # Restructure the ball tree so the points are in locally
        # contiguous blocks of memory (local by branch + leaf).
        if reorder: self._fix_order(self.tree, self.sq_sums, self.radii, order)

    # Restructure the ball tree so the points are in locally
    # contiguous blocks of memory (local by branch + leaf).
    def reorder(self):
        order = self._get_order(search=True)
        self._fix_order(self.tree, self.sq_sums, self.radii, order)

    # Find the "k" nearest neighbors to all points in z. Uses the same
    # interface as the "BallTree.nearest" function, see help for more info.
    def query(self, *args, **kwargs): return self.nearest(*args, **kwargs)

    # Find the "k" nearest neighbors to all points in z.
    def nearest(self, z, k=1, leaf_size=None, return_distance=True,
                transpose=True, max_search=None, look_ahead=None,
                randomized=False):
        # Get the leaf size.
        if (leaf_size is None): leaf_size = self.leaf_size
        # If only a single point was given, convert it to a matrix.
        if (len(z.shape) == 1): z = np.array([z])
        # Transpose the points if appropriate.
        if transpose: z = z.T
        # Make sure the 'k' value isn't bigger than the tree size.
        k = min(k, self.built)
        # Initialize holders for output.
        points  = np.asarray(z, order='F', dtype=self.ttype)
        indices = np.ones((k, points.shape[1]), order='F', dtype='int64')
        dists   = np.ones((k, points.shape[1]), order='F', dtype='float64')
        order = self._get_order(search=True)
        # Compute the nearest neighbors.
        if look_ahead is None:
            self._bt_nearest(points, k, self.tree, self.sq_sums, self.radii,
                             order, leaf_size, indices, dists,
                             to_search=max_search)
        else:
            self._bt_approx_nearest(points, k, self.tree, self.sq_sums, self.radii,
                                    order, leaf_size, indices, dists,
                                    look_ahead=look_ahead, randomized=randomized)
        # Return the results.
        if return_distance:
            if transpose: return dists.T, indices.T - 1
            else:         return dists,   indices - 1
        else:
            if transpose: return indices.T - 1
            else:         return indices - 1

    # Summarize this tree.
    def __str__(self):
        if not hasattr(self,"tree"): return "empty BallTree"
        return f"BallTree {self.tree.shape[::-1]}[:{self.size}] -- {self.built} built"

    # Return the usable length of this tree.
    def __len__(self):
        if not hasattr(self,"built"): return 0
        else:                         return self.built

    # Get an index point from the tree.
    def __getitem__(self, index):
        return self.tree[:,self.order[:self.built][index]-1]

    # Prune this tree and compact its points into the front of the
    # array, adjust '.size' and '.built' accordingly.
    def prune(self, levels=1, full=True, build=True, method="root"):
        assert(levels >= 1)
        assert(len(self) > 0)
        # Build the tree out if it is not fully built.
        if (full and (not (self.built == self.size))): self.build()
        # Get the size of the built portion of the tree.
        size = self.built
        # Compute the indices of the inner children (50th percentile) points.
        if (levels == 1):
            # Use the middle child for 1 level.
            indices = np.array([1], dtype=np.int64)
        else:
            # Handle the different methods of pruning the tree.
            if method in {"inner", "outer"}:
                to_keep = min(size, 1 + 2**(levels-1))
                # Otherwise, get the root, first outer, and all inners.
                indices = np.zeros(to_keep, dtype=np.int64)
                # Get the indices of the inner children of the built tree.
                if method == "inner":
                    indices[0] = 1
                    indices[1] = min((size + 2) // 2, size-1) + 1
                    self._inner(size, levels, indices[2:])
                elif method == "outer":
                    indices[0] = 1
                    indices[1] = 2
                    self._outer(size, levels, indices[2:])
            else:
                to_keep = min(size, 2**levels - 1)
                # Simply grab the root of the tree.
                indices = np.zeros(to_keep, dtype=np.int64)
                self._top(size, levels, indices)
            indices[:] -= 1
        # Stop this operation if the tree will remain unchanged.
        if (len(indices) == size): return
        # Adjust the recorded size and amount built of this tree.
        self.size = len(indices)
        self.built = 1
        # Keep those indices only.
        self.order[:self.size] = self.order[indices]
        # Rebuild the tree if desired.
        if build: self.build(root=0)


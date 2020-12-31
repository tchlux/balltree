'''This Python code is an automatically generated wrapper
for Fortran code made by 'fmodpy'. The original documentation
for the Fortran source code follows.


'''

import os
import ctypes
import numpy

# --------------------------------------------------------------------
#               CONFIGURATION
# 
_verbose = True
_fort_compiler = "gfortran"
_shared_object_name = "ball_tree.so"
_this_directory = os.path.dirname(os.path.abspath(__file__))
_path_to_lib = os.path.join(_this_directory, _shared_object_name)
_compile_options = ['-fPIC', '-shared', '-O3', '-fopenmp']
_ordered_dependencies = ['swap.f90', 'prune.f90', 'fast_select.f90', 'fast_sort.f90', 'ball_tree.f90', 'ball_tree_c_wrapper.f90']
# 
# --------------------------------------------------------------------
#               AUTO-COMPILING
#
# Try to import the existing object. If that fails, recompile and then try.
try:
    clib = ctypes.CDLL(_path_to_lib)
except:
    # Remove the shared object if it exists, because it is faulty.
    if os.path.exists(_shared_object_name):
        os.remove(_shared_object_name)
    # Compile a new shared object.
    _command = " ".join([_fort_compiler] + _compile_options + ["-o", _shared_object_name] + _ordered_dependencies)
    if _verbose:
        print("Running system command with arguments")
        print("  ", _command)
    # Run the compilation command.
    import subprocess
    subprocess.run(_command, shell=True, cwd=_this_directory)
    # Import the shared object file as a C library with ctypes.
    clib = ctypes.CDLL(_path_to_lib)
# --------------------------------------------------------------------


class ball_tree_r64:
    ''''''

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine BUILD_TREE
    
    def build_tree(self, points, sq_sums, radii, order, root=None, leaf_size=None, computed_sq_sums=None):
        '''! Re-arrange elements of POINTS into a binary ball tree.'''
        
        # Setting up "points"
        if ((not issubclass(type(points), numpy.ndarray)) or
            (not numpy.asarray(points).flags.f_contiguous) or
            (not (points.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'points' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            points = numpy.asarray(points, dtype=ctypes.c_double, order='F')
        points_dim_1 = ctypes.c_int(points.shape[0])
        points_dim_2 = ctypes.c_int(points.shape[1])
        
        # Setting up "sq_sums"
        sq_sums_dim_1 = ctypes.c_int(sq_sums.shape[0])
        
        # Setting up "radii"
        radii_dim_1 = ctypes.c_int(radii.shape[0])
        
        # Setting up "order"
        if ((not issubclass(type(order), numpy.ndarray)) or
            (not numpy.asarray(order).flags.f_contiguous) or
            (not (order.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'order' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            order = numpy.asarray(order, dtype=ctypes.c_long, order='F')
        order_dim_1 = ctypes.c_int(order.shape[0])
        
        # Setting up "root"
        root_present = ctypes.c_bool(True)
        if (root is None):
            root_present = ctypes.c_bool(False)
            root = ctypes.c_long()
        if (type(root) is not ctypes.c_long): root = ctypes.c_long(root)
        
        # Setting up "leaf_size"
        leaf_size_present = ctypes.c_bool(True)
        if (leaf_size is None):
            leaf_size_present = ctypes.c_bool(False)
            leaf_size = ctypes.c_long()
        if (type(leaf_size) is not ctypes.c_long): leaf_size = ctypes.c_long(leaf_size)
        
        # Setting up "computed_sq_sums"
        computed_sq_sums_present = ctypes.c_bool(True)
        if (computed_sq_sums is None):
            computed_sq_sums_present = ctypes.c_bool(False)
            computed_sq_sums = ctypes.c_bool()
        if (type(computed_sq_sums) is not ctypes.c_bool): computed_sq_sums = ctypes.c_bool(computed_sq_sums)
    
        # Call C-accessible Fortran wrapper.
        clib.c_build_tree(ctypes.byref(points_dim_1), ctypes.byref(points_dim_2), ctypes.c_void_p(points.ctypes.data), ctypes.byref(sq_sums_dim_1), ctypes.c_void_p(sq_sums.ctypes.data), ctypes.byref(radii_dim_1), ctypes.c_void_p(radii.ctypes.data), ctypes.byref(order_dim_1), ctypes.c_void_p(order.ctypes.data), ctypes.byref(root_present), ctypes.byref(root), ctypes.byref(leaf_size_present), ctypes.byref(leaf_size), ctypes.byref(computed_sq_sums_present), ctypes.byref(computed_sq_sums))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return points, sq_sums, radii, order

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine NEAREST
    
    def nearest(self, points, k, tree, sq_sums, radii, order, leaf_size, indices=None, dists=None, to_search=None):
        '''! Compute the K nearest elements of TREE to each point in POINTS.'''
        
        # Setting up "points"
        if ((not issubclass(type(points), numpy.ndarray)) or
            (not numpy.asarray(points).flags.f_contiguous) or
            (not (points.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'points' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            points = numpy.asarray(points, dtype=ctypes.c_double, order='F')
        points_dim_1 = ctypes.c_int(points.shape[0])
        points_dim_2 = ctypes.c_int(points.shape[1])
        
        # Setting up "k"
        if (type(k) is not ctypes.c_long): k = ctypes.c_long(k)
        
        # Setting up "tree"
        if ((not issubclass(type(tree), numpy.ndarray)) or
            (not numpy.asarray(tree).flags.f_contiguous) or
            (not (tree.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'tree' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            tree = numpy.asarray(tree, dtype=ctypes.c_double, order='F')
        tree_dim_1 = ctypes.c_int(tree.shape[0])
        tree_dim_2 = ctypes.c_int(tree.shape[1])
        
        # Setting up "sq_sums"
        if ((not issubclass(type(sq_sums), numpy.ndarray)) or
            (not numpy.asarray(sq_sums).flags.f_contiguous) or
            (not (sq_sums.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'sq_sums' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            sq_sums = numpy.asarray(sq_sums, dtype=ctypes.c_double, order='F')
        sq_sums_dim_1 = ctypes.c_int(sq_sums.shape[0])
        
        # Setting up "radii"
        if ((not issubclass(type(radii), numpy.ndarray)) or
            (not numpy.asarray(radii).flags.f_contiguous) or
            (not (radii.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'radii' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            radii = numpy.asarray(radii, dtype=ctypes.c_double, order='F')
        radii_dim_1 = ctypes.c_int(radii.shape[0])
        
        # Setting up "order"
        if ((not issubclass(type(order), numpy.ndarray)) or
            (not numpy.asarray(order).flags.f_contiguous) or
            (not (order.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'order' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            order = numpy.asarray(order, dtype=ctypes.c_long, order='F')
        order_dim_1 = ctypes.c_int(order.shape[0])
        
        # Setting up "leaf_size"
        if (type(leaf_size) is not ctypes.c_long): leaf_size = ctypes.c_long(leaf_size)
        
        # Setting up "indices"
        if (indices is None):
            indices = numpy.zeros(shape=(k, points.shape[1]), dtype=ctypes.c_long, order='F')
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        indices_dim_2 = ctypes.c_int(indices.shape[1])
        
        # Setting up "dists"
        if (dists is None):
            dists = numpy.zeros(shape=(k, points.shape[1]), dtype=ctypes.c_double, order='F')
        dists_dim_1 = ctypes.c_int(dists.shape[0])
        dists_dim_2 = ctypes.c_int(dists.shape[1])
        
        # Setting up "to_search"
        to_search_present = ctypes.c_bool(True)
        if (to_search is None):
            to_search_present = ctypes.c_bool(False)
            to_search = ctypes.c_long()
        if (type(to_search) is not ctypes.c_long): to_search = ctypes.c_long(to_search)
    
        # Call C-accessible Fortran wrapper.
        clib.c_nearest(ctypes.byref(points_dim_1), ctypes.byref(points_dim_2), ctypes.c_void_p(points.ctypes.data), ctypes.byref(k), ctypes.byref(tree_dim_1), ctypes.byref(tree_dim_2), ctypes.c_void_p(tree.ctypes.data), ctypes.byref(sq_sums_dim_1), ctypes.c_void_p(sq_sums.ctypes.data), ctypes.byref(radii_dim_1), ctypes.c_void_p(radii.ctypes.data), ctypes.byref(order_dim_1), ctypes.c_void_p(order.ctypes.data), ctypes.byref(leaf_size), ctypes.byref(indices_dim_1), ctypes.byref(indices_dim_2), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(dists_dim_1), ctypes.byref(dists_dim_2), ctypes.c_void_p(dists.ctypes.data), ctypes.byref(to_search_present), ctypes.byref(to_search))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return indices, dists

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine PT_NEAREST
    
    def pt_nearest(self, point, k, tree, sq_sums, radii, order, leaf_size, indices, dists, checks=None, found=None, pt_ss=None):
        '''! Compute the K nearest elements of TREE to each point in POINTS.'''
        
        # Setting up "point"
        if ((not issubclass(type(point), numpy.ndarray)) or
            (not numpy.asarray(point).flags.f_contiguous) or
            (not (point.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'point' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            point = numpy.asarray(point, dtype=ctypes.c_double, order='F')
        point_dim_1 = ctypes.c_int(point.shape[0])
        
        # Setting up "k"
        if (type(k) is not ctypes.c_long): k = ctypes.c_long(k)
        
        # Setting up "tree"
        if ((not issubclass(type(tree), numpy.ndarray)) or
            (not numpy.asarray(tree).flags.f_contiguous) or
            (not (tree.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'tree' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            tree = numpy.asarray(tree, dtype=ctypes.c_double, order='F')
        tree_dim_1 = ctypes.c_int(tree.shape[0])
        tree_dim_2 = ctypes.c_int(tree.shape[1])
        
        # Setting up "sq_sums"
        if ((not issubclass(type(sq_sums), numpy.ndarray)) or
            (not numpy.asarray(sq_sums).flags.f_contiguous) or
            (not (sq_sums.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'sq_sums' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            sq_sums = numpy.asarray(sq_sums, dtype=ctypes.c_double, order='F')
        sq_sums_dim_1 = ctypes.c_int(sq_sums.shape[0])
        
        # Setting up "radii"
        if ((not issubclass(type(radii), numpy.ndarray)) or
            (not numpy.asarray(radii).flags.f_contiguous) or
            (not (radii.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'radii' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            radii = numpy.asarray(radii, dtype=ctypes.c_double, order='F')
        radii_dim_1 = ctypes.c_int(radii.shape[0])
        
        # Setting up "order"
        if ((not issubclass(type(order), numpy.ndarray)) or
            (not numpy.asarray(order).flags.f_contiguous) or
            (not (order.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'order' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            order = numpy.asarray(order, dtype=ctypes.c_long, order='F')
        order_dim_1 = ctypes.c_int(order.shape[0])
        
        # Setting up "leaf_size"
        if (type(leaf_size) is not ctypes.c_long): leaf_size = ctypes.c_long(leaf_size)
        
        # Setting up "indices"
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "dists"
        dists_dim_1 = ctypes.c_int(dists.shape[0])
        
        # Setting up "checks"
        checks_present = ctypes.c_bool(True)
        if (checks is None):
            checks_present = ctypes.c_bool(False)
            checks = ctypes.c_long()
        if (type(checks) is not ctypes.c_long): checks = ctypes.c_long(checks)
        
        # Setting up "found"
        found_present = ctypes.c_bool(True)
        if (found is None):
            found_present = ctypes.c_bool(False)
            found = ctypes.c_long()
        if (type(found) is not ctypes.c_long): found = ctypes.c_long(found)
        
        # Setting up "pt_ss"
        pt_ss_present = ctypes.c_bool(True)
        if (pt_ss is None):
            pt_ss_present = ctypes.c_bool(False)
            pt_ss = ctypes.c_double()
        if (type(pt_ss) is not ctypes.c_double): pt_ss = ctypes.c_double(pt_ss)
    
        # Call C-accessible Fortran wrapper.
        clib.c_pt_nearest(ctypes.byref(point_dim_1), ctypes.c_void_p(point.ctypes.data), ctypes.byref(k), ctypes.byref(tree_dim_1), ctypes.byref(tree_dim_2), ctypes.c_void_p(tree.ctypes.data), ctypes.byref(sq_sums_dim_1), ctypes.c_void_p(sq_sums.ctypes.data), ctypes.byref(radii_dim_1), ctypes.c_void_p(radii.ctypes.data), ctypes.byref(order_dim_1), ctypes.c_void_p(order.ctypes.data), ctypes.byref(leaf_size), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(dists_dim_1), ctypes.c_void_p(dists.ctypes.data), ctypes.byref(checks_present), ctypes.byref(checks), ctypes.byref(found_present), ctypes.byref(found), ctypes.byref(pt_ss_present), ctypes.byref(pt_ss))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return indices, dists, (checks.value if checks_present else None), (found.value if found_present else None)

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine APPROX_NEAREST
    
    def approx_nearest(self, points, k, tree, sq_sums, radii, order, leaf_size, indices=None, dists=None, look_ahead=None):
        '''! Compute the K nearest elements of TREE to each point in POINTS.'''
        
        # Setting up "points"
        if ((not issubclass(type(points), numpy.ndarray)) or
            (not numpy.asarray(points).flags.f_contiguous) or
            (not (points.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'points' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            points = numpy.asarray(points, dtype=ctypes.c_double, order='F')
        points_dim_1 = ctypes.c_int(points.shape[0])
        points_dim_2 = ctypes.c_int(points.shape[1])
        
        # Setting up "k"
        if (type(k) is not ctypes.c_long): k = ctypes.c_long(k)
        
        # Setting up "tree"
        if ((not issubclass(type(tree), numpy.ndarray)) or
            (not numpy.asarray(tree).flags.f_contiguous) or
            (not (tree.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'tree' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            tree = numpy.asarray(tree, dtype=ctypes.c_double, order='F')
        tree_dim_1 = ctypes.c_int(tree.shape[0])
        tree_dim_2 = ctypes.c_int(tree.shape[1])
        
        # Setting up "sq_sums"
        if ((not issubclass(type(sq_sums), numpy.ndarray)) or
            (not numpy.asarray(sq_sums).flags.f_contiguous) or
            (not (sq_sums.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'sq_sums' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            sq_sums = numpy.asarray(sq_sums, dtype=ctypes.c_double, order='F')
        sq_sums_dim_1 = ctypes.c_int(sq_sums.shape[0])
        
        # Setting up "radii"
        if ((not issubclass(type(radii), numpy.ndarray)) or
            (not numpy.asarray(radii).flags.f_contiguous) or
            (not (radii.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'radii' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            radii = numpy.asarray(radii, dtype=ctypes.c_double, order='F')
        radii_dim_1 = ctypes.c_int(radii.shape[0])
        
        # Setting up "order"
        if ((not issubclass(type(order), numpy.ndarray)) or
            (not numpy.asarray(order).flags.f_contiguous) or
            (not (order.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'order' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            order = numpy.asarray(order, dtype=ctypes.c_long, order='F')
        order_dim_1 = ctypes.c_int(order.shape[0])
        
        # Setting up "leaf_size"
        if (type(leaf_size) is not ctypes.c_long): leaf_size = ctypes.c_long(leaf_size)
        
        # Setting up "indices"
        if (indices is None):
            indices = numpy.zeros(shape=(k, points.shape[1]), dtype=ctypes.c_long, order='F')
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        indices_dim_2 = ctypes.c_int(indices.shape[1])
        
        # Setting up "dists"
        if (dists is None):
            dists = numpy.zeros(shape=(k, points.shape[1]), dtype=ctypes.c_double, order='F')
        dists_dim_1 = ctypes.c_int(dists.shape[0])
        dists_dim_2 = ctypes.c_int(dists.shape[1])
        
        # Setting up "look_ahead"
        look_ahead_present = ctypes.c_bool(True)
        if (look_ahead is None):
            look_ahead_present = ctypes.c_bool(False)
            look_ahead = ctypes.c_long()
        if (type(look_ahead) is not ctypes.c_long): look_ahead = ctypes.c_long(look_ahead)
    
        # Call C-accessible Fortran wrapper.
        clib.c_approx_nearest(ctypes.byref(points_dim_1), ctypes.byref(points_dim_2), ctypes.c_void_p(points.ctypes.data), ctypes.byref(k), ctypes.byref(tree_dim_1), ctypes.byref(tree_dim_2), ctypes.c_void_p(tree.ctypes.data), ctypes.byref(sq_sums_dim_1), ctypes.c_void_p(sq_sums.ctypes.data), ctypes.byref(radii_dim_1), ctypes.c_void_p(radii.ctypes.data), ctypes.byref(order_dim_1), ctypes.c_void_p(order.ctypes.data), ctypes.byref(leaf_size), ctypes.byref(indices_dim_1), ctypes.byref(indices_dim_2), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(dists_dim_1), ctypes.byref(dists_dim_2), ctypes.c_void_p(dists.ctypes.data), ctypes.byref(look_ahead_present), ctypes.byref(look_ahead))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return indices, dists

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine PT_APPROX_NEAREST
    
    def pt_approx_nearest(self, point, k, tree, sq_sums, radii, order, leaf_size, indices, dists, look_ahead, found=None, pt_ss=None):
        '''! Compute the K nearest elements of TREE to each point in POINTS
! using the "look ahead" method for determining tree traversal.'''
        
        # Setting up "point"
        if ((not issubclass(type(point), numpy.ndarray)) or
            (not numpy.asarray(point).flags.f_contiguous) or
            (not (point.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'point' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            point = numpy.asarray(point, dtype=ctypes.c_double, order='F')
        point_dim_1 = ctypes.c_int(point.shape[0])
        
        # Setting up "k"
        if (type(k) is not ctypes.c_long): k = ctypes.c_long(k)
        
        # Setting up "tree"
        if ((not issubclass(type(tree), numpy.ndarray)) or
            (not numpy.asarray(tree).flags.f_contiguous) or
            (not (tree.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'tree' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            tree = numpy.asarray(tree, dtype=ctypes.c_double, order='F')
        tree_dim_1 = ctypes.c_int(tree.shape[0])
        tree_dim_2 = ctypes.c_int(tree.shape[1])
        
        # Setting up "sq_sums"
        if ((not issubclass(type(sq_sums), numpy.ndarray)) or
            (not numpy.asarray(sq_sums).flags.f_contiguous) or
            (not (sq_sums.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'sq_sums' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            sq_sums = numpy.asarray(sq_sums, dtype=ctypes.c_double, order='F')
        sq_sums_dim_1 = ctypes.c_int(sq_sums.shape[0])
        
        # Setting up "radii"
        if ((not issubclass(type(radii), numpy.ndarray)) or
            (not numpy.asarray(radii).flags.f_contiguous) or
            (not (radii.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'radii' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            radii = numpy.asarray(radii, dtype=ctypes.c_double, order='F')
        radii_dim_1 = ctypes.c_int(radii.shape[0])
        
        # Setting up "order"
        if ((not issubclass(type(order), numpy.ndarray)) or
            (not numpy.asarray(order).flags.f_contiguous) or
            (not (order.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'order' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            order = numpy.asarray(order, dtype=ctypes.c_long, order='F')
        order_dim_1 = ctypes.c_int(order.shape[0])
        
        # Setting up "leaf_size"
        if (type(leaf_size) is not ctypes.c_long): leaf_size = ctypes.c_long(leaf_size)
        
        # Setting up "indices"
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "dists"
        dists_dim_1 = ctypes.c_int(dists.shape[0])
        
        # Setting up "look_ahead"
        if (type(look_ahead) is not ctypes.c_long): look_ahead = ctypes.c_long(look_ahead)
        
        # Setting up "found"
        found_present = ctypes.c_bool(True)
        if (found is None):
            found_present = ctypes.c_bool(False)
            found = ctypes.c_long()
        if (type(found) is not ctypes.c_long): found = ctypes.c_long(found)
        
        # Setting up "pt_ss"
        pt_ss_present = ctypes.c_bool(True)
        if (pt_ss is None):
            pt_ss_present = ctypes.c_bool(False)
            pt_ss = ctypes.c_double()
        if (type(pt_ss) is not ctypes.c_double): pt_ss = ctypes.c_double(pt_ss)
    
        # Call C-accessible Fortran wrapper.
        clib.c_pt_approx_nearest(ctypes.byref(point_dim_1), ctypes.c_void_p(point.ctypes.data), ctypes.byref(k), ctypes.byref(tree_dim_1), ctypes.byref(tree_dim_2), ctypes.c_void_p(tree.ctypes.data), ctypes.byref(sq_sums_dim_1), ctypes.c_void_p(sq_sums.ctypes.data), ctypes.byref(radii_dim_1), ctypes.c_void_p(radii.ctypes.data), ctypes.byref(order_dim_1), ctypes.c_void_p(order.ctypes.data), ctypes.byref(leaf_size), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(dists_dim_1), ctypes.c_void_p(dists.ctypes.data), ctypes.byref(look_ahead), ctypes.byref(found_present), ctypes.byref(found), ctypes.byref(pt_ss_present), ctypes.byref(pt_ss))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return indices, dists, (found.value if found_present else None)

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine FIX_ORDER
    
    def fix_order(self, points, sq_sums, radii, order, copy=None):
        '''! Re-organize a built tree so that it is more usefully packed in memory.'''
        
        # Setting up "points"
        if ((not issubclass(type(points), numpy.ndarray)) or
            (not numpy.asarray(points).flags.f_contiguous) or
            (not (points.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'points' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            points = numpy.asarray(points, dtype=ctypes.c_double, order='F')
        points_dim_1 = ctypes.c_int(points.shape[0])
        points_dim_2 = ctypes.c_int(points.shape[1])
        
        # Setting up "sq_sums"
        sq_sums_dim_1 = ctypes.c_int(sq_sums.shape[0])
        
        # Setting up "radii"
        radii_dim_1 = ctypes.c_int(radii.shape[0])
        
        # Setting up "order"
        if ((not issubclass(type(order), numpy.ndarray)) or
            (not numpy.asarray(order).flags.f_contiguous) or
            (not (order.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'order' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            order = numpy.asarray(order, dtype=ctypes.c_long, order='F')
        order_dim_1 = ctypes.c_int(order.shape[0])
        
        # Setting up "copy"
        copy_present = ctypes.c_bool(True)
        if (copy is None):
            copy_present = ctypes.c_bool(False)
            copy = ctypes.c_bool()
        if (type(copy) is not ctypes.c_bool): copy = ctypes.c_bool(copy)
    
        # Call C-accessible Fortran wrapper.
        clib.c_fix_order(ctypes.byref(points_dim_1), ctypes.byref(points_dim_2), ctypes.c_void_p(points.ctypes.data), ctypes.byref(sq_sums_dim_1), ctypes.c_void_p(sq_sums.ctypes.data), ctypes.byref(radii_dim_1), ctypes.c_void_p(radii.ctypes.data), ctypes.byref(order_dim_1), ctypes.c_void_p(order.ctypes.data), ctypes.byref(copy_present), ctypes.byref(copy))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return points, sq_sums, radii, order

ball_tree_r64 = ball_tree_r64()


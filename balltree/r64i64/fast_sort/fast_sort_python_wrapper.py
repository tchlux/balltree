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
_shared_object_name = "fast_sort.so"
_this_directory = os.path.dirname(os.path.abspath(__file__))
_path_to_lib = os.path.join(_this_directory, _shared_object_name)
_compile_options = ['-fPIC', '-shared', '-O3', '-fopenmp']
_ordered_dependencies = ['swap.f90', 'fast_sort.f90', 'fast_sort_c_wrapper.f90']
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


class fast_sort:
    ''''''

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine ARGSORT
    
    def argsort(self, values, indices, min_size=None):
        '''! ------------------------------------------------------------------
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
!'''
        
        # Setting up "values"
        if ((not issubclass(type(values), numpy.ndarray)) or
            (not numpy.asarray(values).flags.f_contiguous) or
            (not (values.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'values' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            values = numpy.asarray(values, dtype=ctypes.c_double, order='F')
        values_dim_1 = ctypes.c_int(values.shape[0])
        
        # Setting up "indices"
        if ((not issubclass(type(indices), numpy.ndarray)) or
            (not numpy.asarray(indices).flags.f_contiguous) or
            (not (indices.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'indices' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            indices = numpy.asarray(indices, dtype=ctypes.c_long, order='F')
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "min_size"
        min_size_present = ctypes.c_bool(True)
        if (min_size is None):
            min_size_present = ctypes.c_bool(False)
            min_size = ctypes.c_long()
        if (type(min_size) is not ctypes.c_long): min_size = ctypes.c_long(min_size)
    
        # Call C-accessible Fortran wrapper.
        clib.c_argsort(ctypes.byref(values_dim_1), ctypes.c_void_p(values.ctypes.data), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(min_size_present), ctypes.byref(min_size))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return values, indices

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine INSERTION_ARGSORT
    
    def insertion_argsort(self, values, indices):
        '''! Insertion sort (best for small lists).'''
        
        # Setting up "values"
        if ((not issubclass(type(values), numpy.ndarray)) or
            (not numpy.asarray(values).flags.f_contiguous) or
            (not (values.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'values' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            values = numpy.asarray(values, dtype=ctypes.c_double, order='F')
        values_dim_1 = ctypes.c_int(values.shape[0])
        
        # Setting up "indices"
        if ((not issubclass(type(indices), numpy.ndarray)) or
            (not numpy.asarray(indices).flags.f_contiguous) or
            (not (indices.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'indices' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            indices = numpy.asarray(indices, dtype=ctypes.c_long, order='F')
        indices_dim_1 = ctypes.c_int(indices.shape[0])
    
        # Call C-accessible Fortran wrapper.
        clib.c_insertion_argsort(ctypes.byref(values_dim_1), ctypes.c_void_p(values.ctypes.data), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return values, indices

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine ARGPARTITION
    
    def argpartition(self, values, indices):
        '''! This function efficiently partitions values based on the median
! of the first, middle, and last elements of the VALUES array. This
! function returns the index of the pivot.'''
        
        # Setting up "values"
        if ((not issubclass(type(values), numpy.ndarray)) or
            (not numpy.asarray(values).flags.f_contiguous) or
            (not (values.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'values' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            values = numpy.asarray(values, dtype=ctypes.c_double, order='F')
        values_dim_1 = ctypes.c_int(values.shape[0])
        
        # Setting up "indices"
        if ((not issubclass(type(indices), numpy.ndarray)) or
            (not numpy.asarray(indices).flags.f_contiguous) or
            (not (indices.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'indices' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            indices = numpy.asarray(indices, dtype=ctypes.c_long, order='F')
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "left"
        left = ctypes.c_long()
    
        # Call C-accessible Fortran wrapper.
        clib.c_argpartition(ctypes.byref(values_dim_1), ctypes.c_void_p(values.ctypes.data), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(left))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return values, indices, left.value

fast_sort = fast_sort()


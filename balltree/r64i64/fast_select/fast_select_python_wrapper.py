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
_shared_object_name = "fast_select.so"
_this_directory = os.path.dirname(os.path.abspath(__file__))
_path_to_lib = os.path.join(_this_directory, _shared_object_name)
_compile_options = ['-fPIC', '-shared', '-O3', '-fopenmp']
_ordered_dependencies = ['swap.f90', 'fast_select.f90', 'fast_select_c_wrapper.f90']
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


class fast_select:
    ''''''

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine ARGSELECT_R64
    
    def argselect_r64(self, values, indices, k, divisor=None, max_size=None):
        '''! ------------------------------------------------------------------
!                       FastSelect method
!
! Given VALUES list of numbers, rearrange the elements of VALUES
! such that the element at index K has rank K (holds its same
! location as if all of VALUES were sorted). Symmetrically rearrange
! array INDICES to keep track of prior indices.
!
! This algorithm uses the same conceptual approach as Floyd-Rivest,
! but instead of standard-deviation based selection of bounds for
! recursion, a rank-based method is used to pick the subset of
! values that is searched. This simplifies the code and improves
! interpretability, while achieving the same tunable performance.
!
! Arguments:
!
!   VALUES   --  A 1D array of real numbers.
!   INDICES  --  A 1D array of original indices for elements of VALUES.
!   K        --  A positive integer for the rank index about which
!                VALUES should be rearranged.
! Optional:
!
!   DIVISOR  --  A positive integer >= 2 that represents the
!                division factor used for large VALUES arrays.
!   MAX_SIZE --  An integer >= DIVISOR that represents the largest
!                sized VALUES for which the worst-case pivot value
!                selection is tolerable. A worst-case pivot causes
!                O( SIZE(VALUES)^2 ) runtime. This value should be
!                determined heuristically based on compute hardware.
!
! Output:
!
!   The elements of the array VALUES are rearranged such that the
!   element at position VALUES(K) is in the same location it would
!   be if all of VALUES were in sorted order. Also known as,
!   VALUES(K) has rank K.
!
! Arguments'''
        
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
        
        # Setting up "k"
        if (type(k) is not ctypes.c_long): k = ctypes.c_long(k)
        
        # Setting up "divisor"
        divisor_present = ctypes.c_bool(True)
        if (divisor is None):
            divisor_present = ctypes.c_bool(False)
            divisor = ctypes.c_long()
        if (type(divisor) is not ctypes.c_long): divisor = ctypes.c_long(divisor)
        
        # Setting up "max_size"
        max_size_present = ctypes.c_bool(True)
        if (max_size is None):
            max_size_present = ctypes.c_bool(False)
            max_size = ctypes.c_long()
        if (type(max_size) is not ctypes.c_long): max_size = ctypes.c_long(max_size)
    
        # Call C-accessible Fortran wrapper.
        clib.c_argselect_r64(ctypes.byref(values_dim_1), ctypes.c_void_p(values.ctypes.data), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(k), ctypes.byref(divisor_present), ctypes.byref(divisor), ctypes.byref(max_size_present), ctypes.byref(max_size))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return values, indices

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine ARGSELECT_I64
    
    def argselect_i64(self, values, indices, k, divisor=None, max_size=None):
        '''! Arguments'''
        
        # Setting up "values"
        if ((not issubclass(type(values), numpy.ndarray)) or
            (not numpy.asarray(values).flags.f_contiguous) or
            (not (values.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'values' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            values = numpy.asarray(values, dtype=ctypes.c_long, order='F')
        values_dim_1 = ctypes.c_int(values.shape[0])
        
        # Setting up "indices"
        if ((not issubclass(type(indices), numpy.ndarray)) or
            (not numpy.asarray(indices).flags.f_contiguous) or
            (not (indices.dtype == numpy.dtype(ctypes.c_long)))):
            import warnings
            warnings.warn("The provided argument 'indices' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
            indices = numpy.asarray(indices, dtype=ctypes.c_long, order='F')
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "k"
        if (type(k) is not ctypes.c_long): k = ctypes.c_long(k)
        
        # Setting up "divisor"
        divisor_present = ctypes.c_bool(True)
        if (divisor is None):
            divisor_present = ctypes.c_bool(False)
            divisor = ctypes.c_long()
        if (type(divisor) is not ctypes.c_long): divisor = ctypes.c_long(divisor)
        
        # Setting up "max_size"
        max_size_present = ctypes.c_bool(True)
        if (max_size is None):
            max_size_present = ctypes.c_bool(False)
            max_size = ctypes.c_long()
        if (type(max_size) is not ctypes.c_long): max_size = ctypes.c_long(max_size)
    
        # Call C-accessible Fortran wrapper.
        clib.c_argselect_i64(ctypes.byref(values_dim_1), ctypes.c_void_p(values.ctypes.data), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(k), ctypes.byref(divisor_present), ctypes.byref(divisor), ctypes.byref(max_size_present), ctypes.byref(max_size))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return values, indices

fast_select = fast_select()


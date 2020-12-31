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
_shared_object_name = "prune.so"
_this_directory = os.path.dirname(os.path.abspath(__file__))
_path_to_lib = os.path.join(_this_directory, _shared_object_name)
_compile_options = ['-fPIC', '-shared', '-O3', '-fopenmp']
_ordered_dependencies = ['prune.f90', 'prune_c_wrapper.f90']
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


class prune:
    ''''''

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine INNER
    
    def inner(self, tree_size, levels, indices, current_level=None, starting_index=None, next_index=None):
        '''! Get the indices of all the inner children in this tree (the 50th
! percentile points) up to level 'LEVELS'.'''
        
        # Setting up "tree_size"
        if (type(tree_size) is not ctypes.c_long): tree_size = ctypes.c_long(tree_size)
        
        # Setting up "levels"
        if (type(levels) is not ctypes.c_long): levels = ctypes.c_long(levels)
        
        # Setting up "indices"
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "current_level"
        current_level_present = ctypes.c_bool(True)
        if (current_level is None):
            current_level_present = ctypes.c_bool(False)
            current_level = ctypes.c_long()
        if (type(current_level) is not ctypes.c_long): current_level = ctypes.c_long(current_level)
        
        # Setting up "starting_index"
        starting_index_present = ctypes.c_bool(True)
        if (starting_index is None):
            starting_index_present = ctypes.c_bool(False)
            starting_index = ctypes.c_long()
        if (type(starting_index) is not ctypes.c_long): starting_index = ctypes.c_long(starting_index)
        
        # Setting up "next_index"
        next_index_present = ctypes.c_bool(True)
        if (next_index is None):
            next_index_present = ctypes.c_bool(False)
            next_index = ctypes.c_long()
        if (type(next_index) is not ctypes.c_long): next_index = ctypes.c_long(next_index)
    
        # Call C-accessible Fortran wrapper.
        clib.c_inner(ctypes.byref(tree_size), ctypes.byref(levels), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(current_level_present), ctypes.byref(current_level), ctypes.byref(starting_index_present), ctypes.byref(starting_index), ctypes.byref(next_index_present), ctypes.byref(next_index))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return indices, (next_index.value if next_index_present else None)

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine OUTER
    
    def outer(self, tree_size, levels, indices, current_level=None, starting_index=None, next_index=None):
        '''! Get the indices of all the outer children in this tree (extrema
! points) up to level 'LEVELS'.'''
        
        # Setting up "tree_size"
        if (type(tree_size) is not ctypes.c_long): tree_size = ctypes.c_long(tree_size)
        
        # Setting up "levels"
        if (type(levels) is not ctypes.c_long): levels = ctypes.c_long(levels)
        
        # Setting up "indices"
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "current_level"
        current_level_present = ctypes.c_bool(True)
        if (current_level is None):
            current_level_present = ctypes.c_bool(False)
            current_level = ctypes.c_long()
        if (type(current_level) is not ctypes.c_long): current_level = ctypes.c_long(current_level)
        
        # Setting up "starting_index"
        starting_index_present = ctypes.c_bool(True)
        if (starting_index is None):
            starting_index_present = ctypes.c_bool(False)
            starting_index = ctypes.c_long()
        if (type(starting_index) is not ctypes.c_long): starting_index = ctypes.c_long(starting_index)
        
        # Setting up "next_index"
        next_index_present = ctypes.c_bool(True)
        if (next_index is None):
            next_index_present = ctypes.c_bool(False)
            next_index = ctypes.c_long()
        if (type(next_index) is not ctypes.c_long): next_index = ctypes.c_long(next_index)
    
        # Call C-accessible Fortran wrapper.
        clib.c_outer(ctypes.byref(tree_size), ctypes.byref(levels), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(current_level_present), ctypes.byref(current_level), ctypes.byref(starting_index_present), ctypes.byref(starting_index), ctypes.byref(next_index_present), ctypes.byref(next_index))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return indices, (next_index.value if next_index_present else None)

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine TOP
    
    def top(self, tree_size, levels, indices, current_level=None, starting_index=None, next_index=None):
        '''! Get the indices of all the top children in this tree up to level 'LEVELS'.'''
        
        # Setting up "tree_size"
        if (type(tree_size) is not ctypes.c_long): tree_size = ctypes.c_long(tree_size)
        
        # Setting up "levels"
        if (type(levels) is not ctypes.c_long): levels = ctypes.c_long(levels)
        
        # Setting up "indices"
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "current_level"
        current_level_present = ctypes.c_bool(True)
        if (current_level is None):
            current_level_present = ctypes.c_bool(False)
            current_level = ctypes.c_long()
        if (type(current_level) is not ctypes.c_long): current_level = ctypes.c_long(current_level)
        
        # Setting up "starting_index"
        starting_index_present = ctypes.c_bool(True)
        if (starting_index is None):
            starting_index_present = ctypes.c_bool(False)
            starting_index = ctypes.c_long()
        if (type(starting_index) is not ctypes.c_long): starting_index = ctypes.c_long(starting_index)
        
        # Setting up "next_index"
        next_index_present = ctypes.c_bool(True)
        if (next_index is None):
            next_index_present = ctypes.c_bool(False)
            next_index = ctypes.c_long()
        if (type(next_index) is not ctypes.c_long): next_index = ctypes.c_long(next_index)
    
        # Call C-accessible Fortran wrapper.
        clib.c_top(ctypes.byref(tree_size), ctypes.byref(levels), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(current_level_present), ctypes.byref(current_level), ctypes.byref(starting_index_present), ctypes.byref(starting_index), ctypes.byref(next_index_present), ctypes.byref(next_index))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return indices, (next_index.value if next_index_present else None)

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine LEVEL
    
    def level(self, tree_size, layer, indices, current_level=None, starting_index=None, current_found=None):
        '''! Get the indices of the level "LAYER" of this tree, with 0 being
! the root. List of indices is returned in "INDICES".'''
        
        # Setting up "tree_size"
        if (type(tree_size) is not ctypes.c_long): tree_size = ctypes.c_long(tree_size)
        
        # Setting up "layer"
        if (type(layer) is not ctypes.c_long): layer = ctypes.c_long(layer)
        
        # Setting up "indices"
        indices_dim_1 = ctypes.c_int(indices.shape[0])
        
        # Setting up "found"
        found = ctypes.c_long()
        
        # Setting up "current_level"
        current_level_present = ctypes.c_bool(True)
        if (current_level is None):
            current_level_present = ctypes.c_bool(False)
            current_level = ctypes.c_long()
        if (type(current_level) is not ctypes.c_long): current_level = ctypes.c_long(current_level)
        
        # Setting up "starting_index"
        starting_index_present = ctypes.c_bool(True)
        if (starting_index is None):
            starting_index_present = ctypes.c_bool(False)
            starting_index = ctypes.c_long()
        if (type(starting_index) is not ctypes.c_long): starting_index = ctypes.c_long(starting_index)
        
        # Setting up "current_found"
        current_found_present = ctypes.c_bool(True)
        if (current_found is None):
            current_found_present = ctypes.c_bool(False)
            current_found = ctypes.c_long()
        if (type(current_found) is not ctypes.c_long): current_found = ctypes.c_long(current_found)
    
        # Call C-accessible Fortran wrapper.
        clib.c_level(ctypes.byref(tree_size), ctypes.byref(layer), ctypes.byref(indices_dim_1), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(found), ctypes.byref(current_level_present), ctypes.byref(current_level), ctypes.byref(starting_index_present), ctypes.byref(starting_index), ctypes.byref(current_found_present), ctypes.byref(current_found))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return indices, found.value

prune = prune()


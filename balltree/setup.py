#!python 

# Any python code that should be executed (as if __main__) during
# setup should go here to prepare the module on a new computer.

print()
print("RUNNING SETUP SCRIPT NOW!")
print()

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



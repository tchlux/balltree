#!/bin/zsh
gfortran -c -O3 swap.f90 prune.f90 fast_select.f90 fast_sort.f90 ball_tree.f90
rm *.o

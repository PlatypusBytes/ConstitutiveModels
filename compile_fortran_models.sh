#!/bin/bash

# Create build directories if they don't exist
mkdir -p build_Fortran/lib

# Compile UMAT_LinearElasticity.f90 - first create object file, then shared library
gfortran -c -o build_Fortran/UMAT_LinearElasticity.o fortran_models/linear_elastic/UMAT_LinearElasticity.f90
gfortran -shared -o build_Fortran/lib/linear_elastic.so build_Fortran/UMAT_LinearElasticity.o

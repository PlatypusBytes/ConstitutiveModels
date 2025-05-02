REM compile all fortran models with gfortran, uncomment the line below to compile the models locally

if not exist "build_Fortran" (
    mkdir "build_Fortran"
)
if not exist "build_Fortran\lib" (
    mkdir "build_Fortran\lib"
)

gfortran -c -o build_Fortran/UMAT_LinearElasticity.o fortran_models/linear_elastic/UMAT_LinearElasticity.f90
gfortran -shared -o build_Fortran/lib/linear_elastic.dll build_Fortran/UMAT_LinearElasticity.o



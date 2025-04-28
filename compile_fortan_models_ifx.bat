REM this script should be used in: Intel oneAPI command prompt for Intel 64 for Visual Studio 2022

if not exist "build_Fortran" (
    mkdir "build_Fortran"
)
if not exist "build_Fortran\lib" (
    mkdir "build_Fortran\lib"
)
ifx -c -o build_Fortran\UMAT_LinearElasticity.obj fortran_models\linear_elastic\UMAT_LinearElasticity.f90
ifx -dll -o build_Fortran\lib\linear_elastic.dll build_Fortran\UMAT_LinearElasticity.obj

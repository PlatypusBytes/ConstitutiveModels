REM this script should be used in: Intel oneAPI command prompt for Intel 64 for Visual Studio 2022
ifx /dll /libs:static /threads UMAT_LinearElasticity.f90 /link /out:linear_elastic.dll
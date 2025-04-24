REM this script should be used in: Native Tools Command Prompt for VS 2022
REM for x64 architecture, use x64 Native Tools Command Prompt for VS 2022
@REM call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64

cl /c /EHsc /W4 /O2 linear_elastic_with_vertical_tens_cutoff.c ../utils.c ../globals.c ../elastic_laws/hookes_law.c
link /DLL /OUT:linear_elastic_with_vertical_tens_cutoff.dll utils.obj linear_elastic_with_vertical_tens_cutoff.obj globals.obj hookes_law.obj /DEF:umat.def
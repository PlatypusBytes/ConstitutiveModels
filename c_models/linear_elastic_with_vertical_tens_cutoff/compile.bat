REM this script should be used in: Native Tools Command Prompt for VS 2022
REM for x64 architecture, use x64 Native Tools Command Prompt for VS 2022

cl /c /EHsc /W3 linear_elastic_with_vertical_tens_cutoff.c /Fo:linear_elastic_with_vertical_tens_cutoff.obj
link /DLL /OUT:linear_elastic_with_vertical_tens_cutoff.dll linear_elastic_with_vertical_tens_cutoff.obj /DEF:umat.def
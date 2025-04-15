REM this script should be used in: Native Tools Command Prompt for VS 2022
REM for x64 architecture, use x64 Native Tools Command Prompt for VS 2022

@rem locally this script can also be used in another command prompt when setting up the environment
@REM call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64

@rem debug:
@REM cl /c /EHsc /W3 /Zi /MD matsuoka_nakai.c ../utils.c ../stress_utils.c ../globals.c ../yield_surfaces/matsuoka_nakai_surface.c ../elastic_laws/hookes_law.c

@rem release:
cl /c /EHsc /W3 matsuoka_nakai.c ../utils.c ../stress_utils.c ../globals.c ../yield_surfaces/matsuoka_nakai_surface.c ../elastic_laws/hookes_law.c

@rem link the dll
link /DLL /OUT:matsuoka_nakai.dll matsuoka_nakai.obj utils.obj stress_utils.obj globals.obj matsuoka_nakai_surface.obj hookes_law.obj  /DEF:umat.def
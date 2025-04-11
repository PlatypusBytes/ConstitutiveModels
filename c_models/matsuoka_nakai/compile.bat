REM this script should be used in: Native Tools Command Prompt for VS 2022
REM for x64 architecture, use x64 Native Tools Command Prompt for VS 2022

call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
cl /c /EHsc /W3 matsuoka_nakai.c /Fo:matsuoka_nakai.obj
link /DLL /OUT:matsuoka_nakai.dll matsuoka_nakai.obj /DEF:umat.def
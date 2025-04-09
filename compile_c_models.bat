REM This script compiles all the models in the c_models directory,
@REM REM uncomment the line below to compile the models locally
@REM call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64

cd c_models/linear_elastic_with_vertical_tens_cutoff
call compile.bat
cd ../..

# Constitutive Models
Library of soil constitutive models in UMAT format.

Currently implemented models:

With Fortran ordering (Column major):
- Linear Elastic 
- Mohr-Coulomb with tension cutoff (Hans Teunissen implementation)

With C ordering (Row major):
- Linear Elastic with tension cutoff in vertical direction 
- Linear Elastic with tension cutoff in normal direction for 2D interface elements
- Linear Elastic with tension cutoff in normal direction for 3D interface elements
- Matsuoka-Nakai 


# How to compile
## Fortran models
The models which are implemented in Fortran can be compiled using the following command:

In windows with gfortran:
```bash
    compile_fortran_models.bat
```

In windows with Intel Fortran:
- open 'Intel oneAPI command prompt for Intel 64 for Visual Studio 2022'
- run the following command:
```bash
    compile_fortran_models_ifx.bat
```

In Linux with gfortran:
- run the following command:
```bash
    ./compile_fortran_models.sh

```

## C models
The models which are implemented in C can be compiled using the following command:

In windows with mscv (Microsoft Visual Studio is required):
- open 'x64 Native Tools Command Prompt for VS 2022'
- or open cmd and run the following command: 
```bash
    call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
```
- run the following command:
```bash
    nmake -f Makefile.windows
```

In linux with gcc:
- run the following command:
```bash
    make -f Makefile.linux
```

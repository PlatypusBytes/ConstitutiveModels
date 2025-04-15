
#ifndef GLOBALS_H
#define GLOBALS_H

extern const double PI;        // Define PI as a constant
extern const double ZERO_TOL;  // Define tolerance for zero checks

// use #define, so it can be used in preprocessor directives
#define VOIGTSIZE_3D 6  // Define size of Voigt vector for 3D tensors

// Define tensor component mapping (Voigt, 0-based)
extern const int XX;  // X component
extern const int YY;  // Y component
extern const int ZZ;  // Z component
extern const int XY;  // XY shear component
extern const int YZ;  // YZ shear component
extern const int XZ;  // XZ shear component

#endif  // GLOBALS_H

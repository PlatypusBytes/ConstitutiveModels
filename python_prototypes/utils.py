"""
Python counterpart of c_models/globals.c. The linear algebra of c_models/utils.c is done with numpy.
"""

ZERO_TOL = 1.0e-12
SMALL_VALUE = 1.0e-16

# Voigt ordering (see globals.h): [xx, yy, zz, xy, yz, xz], engineering shear strains
VOIGTSIZE_3D = 6
XX, YY, ZZ, XY, YZ, XZ = range(VOIGTSIZE_3D)

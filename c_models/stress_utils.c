#include <math.h>

#include "globals.h"
#include "stress_utils.h"
#include "utils.h"



void calculate_stress_invariants_3d(const double stress[VOIGTSIZE_3D], double* p, double* J,
                                    double* theta, double* j2, double* j3,
                                    double s_dev[VOIGTSIZE_3D])
{
    // Mean stress pressure p = trace(sigma)/3
    *p = (stress[XX] + stress[YY] + stress[ZZ]) / 3.0;

    // Deviatoric stress tensor s = sigma - p * I
    s_dev[XX] = stress[XX] - *p;
    s_dev[YY] = stress[YY] - *p;
    s_dev[ZZ] = stress[ZZ] - *p;
    s_dev[XY] = stress[XY];
    s_dev[YZ] = stress[YZ];
    s_dev[XZ] = stress[XZ];

    // Second invariant of deviatoric stress J2 = 0.5 * s:s
    // J2 = 0.5 * (sxx^2 + syy^2 + szz^2) + sxy^2 + syz^2 + sxz^2
    *j2 = 0.5 * (s_dev[XX] * s_dev[XX] + s_dev[YY] * s_dev[YY] + s_dev[ZZ] * s_dev[ZZ]) +
          (s_dev[XY] * s_dev[XY] + s_dev[YZ] * s_dev[YZ] + s_dev[XZ] * s_dev[XZ]);

    // for convenience, store sqrt(J2) in J
    *J = sqrt((*j2));

    // Third invariant of deviatoric stress J3 = det(s)
    *j3 = calculate_determinant_voigt_vector_3d(s_dev);

    // Lode angle theta
    if (*J > ZERO_TOL)
    {  // Avoid division by zero if J is very small
        double arg = (sqrt(27.0) / 2.0) * (*j3) / pow(*j2, 1.5);

        // Clamp argument to asin range [-1, 1] due to potential numerical inaccuracies
        double zero_tolerance = 0.0;
        if (arg > 1.0)
        {
            arg = 1.0;
            zero_tolerance = -ZERO_TOL;
        }
        if (arg < -1.0)
        {
            arg = -1.0;
            zero_tolerance = ZERO_TOL;
        }

        *theta = asin(arg) / 3.0 + zero_tolerance;  // Result is in [-pi/6, pi/6]
    }
    else
    {
        *theta = 0.0;  // Set to convention if J is zero
    }
}

void calculate_stress_invariants_derivatives_3d(const double J, const double s_dev[VOIGTSIZE_3D],
                                                const double j2, const double j3,
                                                double dp_dsig[VOIGTSIZE_3D],
                                                double dJ_dsig[VOIGTSIZE_3D],
                                                double dtheta_dsig[VOIGTSIZE_3D])
{
    // Derivatives of invariants w.r.t. sigma (in Voigt)
    // dp/dsigma = [1/3, 1/3, 1/3, 0, 0, 0]^T
    dp_dsig[XX] = 1.0 / 3.0;
    dp_dsig[YY] = 1.0 / 3.0;
    dp_dsig[ZZ] = 1.0 / 3.0;
    dp_dsig[XY] = 0.0;
    dp_dsig[YZ] = 0.0;
    dp_dsig[XZ] = 0.0;

    // calculate dJ/dsigma
    if (J < ZERO_TOL)
    {
        // Handle case where q is very small or zero
        for (int i = 0; i < VOIGTSIZE_3D; ++i)
        {
            dJ_dsig[i] = 0.0;  // Set to zero
        }
    }
    else
    {
        dJ_dsig[XX] = s_dev[XX] / (2.0 * J);
        dJ_dsig[YY] = s_dev[YY] / (2.0 * J);
        dJ_dsig[ZZ] = s_dev[ZZ] / (2.0 * J);
        dJ_dsig[XY] = s_dev[XY] / J;
        dJ_dsig[YZ] = s_dev[YZ] / J;
        dJ_dsig[XZ] = s_dev[XZ] / J;
    }

    // dtheta/dsigma = (dtheta/dJ2)*(dJ2/dsig) + (dtheta/dJ3)*(dJ3/dsig)
    // This is the complex part, requires dJ2/dsig, dJ3/dsig, dtheta/dJ2, dtheta/dJ3

    double dJ2_dsig[VOIGTSIZE_3D] = {s_dev[XX],       s_dev[YY],       s_dev[ZZ],
                                     2.0 * s_dev[XY], 2.0 * s_dev[YZ], 2.0 * s_dev[XZ]};

    double dJ3_dsig[VOIGTSIZE_3D];
    dJ3_dsig[XX] = s_dev[XX] * s_dev[XX] + s_dev[XY] * s_dev[XY] + s_dev[XZ] * s_dev[XZ] -
                   2.0 * j2 / 3.0;  // s_xx^2 + s_xy^2 + s_xz^2 - 2/3 * J2
    dJ3_dsig[YY] = s_dev[YY] * s_dev[YY] + s_dev[XY] * s_dev[XY] + s_dev[YZ] * s_dev[YZ] -
                   2.0 * j2 / 3.0;  // s_yy^2 + s_xy^2 + s_yz^2 - 2/3 * J2
    dJ3_dsig[ZZ] = s_dev[ZZ] * s_dev[ZZ] + s_dev[YZ] * s_dev[YZ] + s_dev[XZ] * s_dev[XZ] -
                   2.0 * j2 / 3.0;  // s_zz^2 + s_yz^2 + s_xz^2 - 2/3 * J2
    dJ3_dsig[XY] = 2.0 * ((s_dev[XX] + s_dev[YY]) * s_dev[XY] +
                          s_dev[XZ] * s_dev[YZ]);  // 2*((s_xx+s_yy)*s_xy + s_xz*s_yz)
    dJ3_dsig[YZ] = 2.0 * ((s_dev[YY] + s_dev[ZZ]) * s_dev[YZ] +
                          s_dev[XY] * s_dev[XZ]);  // 2*((s_yy+s_zz)*s_yz + s_xy*s_xz)
    dJ3_dsig[XZ] = 2.0 * ((s_dev[ZZ] + s_dev[XX]) * s_dev[XZ] +
                          s_dev[XY] * s_dev[YZ]);  // 2*((s_zz+s_xx)*s_xz + s_xy*s_yz)

    double dtheta_dJ2 = 0.0;
    double dtheta_dJ3 = 0.0;

    const double inner_term = 1.0 - (27.0 * j3 * j3) / (4.0 * pow(j2, 3.0));

    // prevent sqrt zero
    if (inner_term > ZERO_TOL)
    {
        dtheta_dJ2 = -pow(3.0, 1.5) * j3 / (4.0 * pow(j2, 2.5) * sqrt(inner_term));
        dtheta_dJ3 = sqrt(3.0) / sqrt(4.0 * pow(j2, 3.0) - 27.0 * j3 * j3);
    }

    // dtheta/dsigma = dtheta_dJ2 * dJ2_dsig + dtheta_dJ3 * dJ3_dsig
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        dtheta_dsig[i] = dtheta_dJ2 * dJ2_dsig[i] + dtheta_dJ3 * dJ3_dsig[i];
    }
}


void calculate_principle_stresses_3d(const double stress[VOIGTSIZE_3D], double principle_stresses[3])
{

    double mean_stress;
    double J;
    double theta;
    double j2;
    double j3;
    double deviatoric_stress[VOIGTSIZE_3D];

    calculate_stress_invariants_3d(stress, &mean_stress,&J, &theta, &j2,&j3, deviatoric_stress);

    // theta is the sine-based Lode angle in [-pi/6, pi/6] (theta = pi/6 in triaxial compression),
    // the cosine roots below require the cosine-based angle in [0, pi/3]
    const double alpha = PI / 6.0 - theta;
    const double radius = 2.0 * J / sqrt(3.0);

    principle_stresses[0] = mean_stress + radius * cos(alpha);
    principle_stresses[1] = mean_stress + radius * cos(alpha - 2.0 * PI / 3.0);
    principle_stresses[2] = mean_stress + radius * cos(alpha + 2.0 * PI / 3.0);

}

void calculate_dq_dsigma_triaxial_state_3d(double dqdsigma[VOIGTSIZE_3D])
{
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        dqdsigma[i] = 0.0;
    }

    dqdsigma[0] = 1;
    dqdsigma[2] = -1;
}




/* ------------------------------------------------------------
   Dot product
   ------------------------------------------------------------ */
static double dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] +
           a[1]*b[1] +
           a[2]*b[2];
}


/* ------------------------------------------------------------
   Norm
   ------------------------------------------------------------ */
static double norm3(const double a[3])
{
    return sqrt(dot3(a, a));
}


/* ------------------------------------------------------------
   Normalize
   ------------------------------------------------------------ */
static void normalize3(double a[3])
{
    double n = norm3(a);

    if (n > 1.0e-15) {
        a[0] /= n;
        a[1] /= n;
        a[2] /= n;
    }
}


/* ------------------------------------------------------------
   Cross product
   ------------------------------------------------------------ */
static void cross3(
    const double a[3],
    const double b[3],
    double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}


/* ------------------------------------------------------------
   Swap eigenvectors/eigenvalues
   ------------------------------------------------------------ */
static void swap_eigenpairs(
    double eig[3],
    double Q[3][3],
    int a,
    int b)
{
    double tmp;

    tmp = eig[a];
    eig[a] = eig[b];
    eig[b] = tmp;

    for (int i = 0; i < 3; i++) {
        tmp = Q[i][a];
        Q[i][a] = Q[i][b];
        Q[i][b] = tmp;
    }
}


/* ------------------------------------------------------------
   Jacobi eigensolver for symmetric 3x3 matrix.

   Input:
       A = symmetric matrix

   Output:
       eig = eigenvalues
       Q   = eigenvectors stored as columns

   eig[0] corresponds to Q[:,0], etc.
   ------------------------------------------------------------ */
static void jacobi_eigen_3x3(
    double A[3][3],
    double eig[3],
    double Q[3][3])
{
    /* Initialize Q = identity */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int iter = 0; iter < JACOBI_MAX_ITER; iter++) {

        /*
         * Find largest off-diagonal element.
         */
        int p = 0;
        int q = 1;

        double max_off = fabs(A[0][1]);

        if (fabs(A[0][2]) > max_off) {
            max_off = fabs(A[0][2]);
            p = 0;
            q = 2;
        }

        if (fabs(A[1][2]) > max_off) {
            max_off = fabs(A[1][2]);
            p = 1;
            q = 2;
        }

        if (max_off < EIG_TOL)
            break;


        /*
         * Jacobi rotation.
         */
        double app = A[p][p];
        double aqq = A[q][q];
        double apq = A[p][q];

        double phi = 0.5 * atan2(
            2.0 * apq,
            aqq - app
        );

        double c = cos(phi);
        double s = sin(phi);


        /*
         * Rotate matrix.
         */
        A[p][p] =
            c*c*app
            - 2.0*c*s*apq
            + s*s*aqq;

        A[q][q] =
            s*s*app
            + 2.0*c*s*apq
            + c*c*aqq;

        A[p][q] = 0.0;
        A[q][p] = 0.0;


        for (int k = 0; k < 3; k++) {

            if (k == p || k == q)
                continue;

            double akp = A[k][p];
            double akq = A[k][q];

            A[k][p] = c*akp - s*akq;
            A[p][k] = A[k][p];

            A[k][q] = s*akp + c*akq;
            A[q][k] = A[k][q];
        }


        /*
         * Update eigenvectors.
         */
        for (int k = 0; k < 3; k++) {

            double qkp = Q[k][p];
            double qkq = Q[k][q];

            Q[k][p] = c*qkp - s*qkq;
            Q[k][q] = s*qkp + c*qkq;
        }
    }


    /*
     * Extract eigenvalues.
     */
    eig[0] = A[0][0];
    eig[1] = A[1][1];
    eig[2] = A[2][2];


    /*
     * Sort descending:
     *
     * sigma1 >= sigma2 >= sigma3
     */
    if (eig[0] < eig[1])
        swap_eigenpairs(eig, Q, 0, 1);

    if (eig[1] < eig[2])
        swap_eigenpairs(eig, Q, 1, 2);

    if (eig[0] < eig[1])
        swap_eigenpairs(eig, Q, 0, 1);
}


static void stress_voigt_matrix(
    const double stress[6],
    double A[3][3])
{
    A[0][0] = stress[0];  // sxx
    A[1][1] = stress[1];  // syy
    A[2][2] = stress[2];  // szz

    A[0][1] = stress[3];  // txy
    A[1][0] = stress[3];

    A[1][2] = stress[4];  // tyz
    A[2][1] = stress[4];

    A[0][2] = stress[5];  // txz
    A[2][0] = stress[5];
}

void calculate_principal_system(
    const double stress[6],
    double principal_stress[3],
    double Q[3][3])
{
    double A[3][3];

    stress_voigt_matrix(stress, A);

    jacobi_eigen_3x3(
        A,
        principal_stress,
        Q
    );
}

void calculate_stress_from_principal_system(
    const double principal_stress[3],
    double Q[3][3],
    double stress[VOIGTSIZE_3D])
{
    for (int i = 0; i < VOIGTSIZE_3D; i++)
        stress[i] = 0.0;

    // stress_ij = sum_k principal_stress_k * Q_ik * Q_jk
    for (int k = 0; k < 3; k++) {
        stress[XX] += principal_stress[k] * Q[0][k] * Q[0][k];
        stress[YY] += principal_stress[k] * Q[1][k] * Q[1][k];
        stress[ZZ] += principal_stress[k] * Q[2][k] * Q[2][k];
        stress[XY] += principal_stress[k] * Q[0][k] * Q[1][k];
        stress[YZ] += principal_stress[k] * Q[1][k] * Q[2][k];
        stress[XZ] += principal_stress[k] * Q[0][k] * Q[2][k];
    }
}
#include <math.h>
#include <stddef.h>

#include "../globals.h"
#include "matsuoka_nakai_surface.h"

/**
 * @brief Lode angle function C(xi) = cos(acos(beta xi) / 3 - gamma pi / 6) and its first and
 * second derivatives with respect to xi = sin(3 theta).
 */
static void calculate_lode_function(const double xi, const MatsuokaNakaiConstants constants,
                                    double* C, double* dC_dxi, double* d2C_dxi2)
{
    const double beta = constants.beta;
    const double w = acos(beta * xi) / 3.0 - constants.gamma * PI / 6.0;
    const double root = sqrt(1.0 - beta * beta * xi * xi);  // > 0 for |beta| < 1

    *C = cos(w);
    *dC_dxi = sin(w) * beta / (3.0 * root);
    *d2C_dxi2 = -cos(w) * beta * beta / (9.0 * root * root) +
                sin(w) * beta * beta * beta * xi / (3.0 * root * root * root);
}

void calculate_yield_function(const double p, const double theta, const double J,
                              const MatsuokaNakaiConstants constants, double* f)
{
    *f = constants.M * p - constants.K +
         J * constants.alpha *
             cos(acos(constants.beta * sin(3.0 * theta)) / 3.0 - constants.gamma * PI / 6.0);

    // in paper (compression positive p):
    //*f  = -(K + M * p) + J * alpha * cos(acos(beta * sin(3 * theta)) / 3 - gamma * PI / 6);
}

void calculate_yield_derivatives(const double s_dev[VOIGTSIZE_3D], const double j2,
                                 const double j3, const MatsuokaNakaiConstants constants,
                                 double grad[VOIGTSIZE_3D],
                                 double hessian[VOIGTSIZE_3D * VOIGTSIZE_3D])
{
    // f = M p - K + h(J2, J3), with h = alpha sqrt(J2) C(xi) and xi = 3 sqrt(3) / 2 J3 / J2^(3/2)
    // grad = M dp/dsig + h_2 dJ2/dsig + h_3 dJ3/dsig
    // hessian = h_2 d2J2 + h_3 d2J3 + h_22 dJ2 dJ2^T + h_23 (dJ2 dJ3^T + dJ3 dJ2^T) + h_33 dJ3 dJ3^T

    const double sxx = s_dev[XX], syy = s_dev[YY], szz = s_dev[ZZ];
    const double sxy = s_dev[XY], syz = s_dev[YZ], sxz = s_dev[XZ];

    for (int i = 0; i < VOIGTSIZE_3D; ++i) grad[i] = (i < 3) ? constants.M / 3.0 : 0.0;
    if (hessian != NULL)
        for (int i = 0; i < VOIGTSIZE_3D * VOIGTSIZE_3D; ++i) hessian[i] = 0.0;

    // on the hydrostatic axis the deviatoric part is not differentiable (apex of the cone)
    if (j2 <= ZERO_TOL * ZERO_TOL) return;

    // dJ2/dsig and dJ3/dsig (shear components are derivatives with respect to the Voigt
    // component, i.e. twice the tensor derivative)
    const double dJ2[VOIGTSIZE_3D] = {sxx, syy, szz, 2.0 * sxy, 2.0 * syz, 2.0 * sxz};
    double dJ3[VOIGTSIZE_3D];
    dJ3[XX] = syy * szz - syz * syz + j2 / 3.0;
    dJ3[YY] = sxx * szz - sxz * sxz + j2 / 3.0;
    dJ3[ZZ] = sxx * syy - sxy * sxy + j2 / 3.0;
    dJ3[XY] = 2.0 * (syz * sxz - szz * sxy);
    dJ3[YZ] = 2.0 * (sxy * sxz - sxx * syz);
    dJ3[XZ] = 2.0 * (sxy * syz - syy * sxz);

    // Lode angle function of xi = sin(3 theta), and derivatives of xi with respect to J2 and J3
    const double J = sqrt(j2);
    const double c0 = 1.5 * sqrt(3.0);
    double xi = c0 * j3 / (j2 * J);
    if (xi > 1.0) xi = 1.0;
    if (xi < -1.0) xi = -1.0;

    double C, C1, C2;
    calculate_lode_function(xi, constants, &C, &C1, &C2);

    const double xi_2 = -1.5 * xi / j2;
    const double xi_3 = c0 / (j2 * J);
    const double xi_22 = 3.75 * xi / (j2 * j2);
    const double xi_23 = -1.5 * xi_3 / j2;

    const double alpha = constants.alpha;
    const double h_2 = alpha * (0.5 * C / J + J * C1 * xi_2);
    const double h_3 = alpha * J * C1 * xi_3;

    for (int i = 0; i < VOIGTSIZE_3D; ++i) grad[i] += h_2 * dJ2[i] + h_3 * dJ3[i];

    if (hessian == NULL) return;

    const double h_22 =
        alpha * (-0.25 * C / (J * j2) + C1 * xi_2 / J + J * (C2 * xi_2 * xi_2 + C1 * xi_22));
    const double h_23 = alpha * (0.5 * C1 * xi_3 / J + J * (C2 * xi_2 * xi_3 + C1 * xi_23));
    const double h_33 = alpha * J * C2 * xi_3 * xi_3;

    // deviatoric projection of the normal components
    double P[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) P[i][j] = (i == j ? 1.0 : 0.0) - 1.0 / 3.0;

    // second derivative of J3 = det(s): normal-normal block P A P, normal-shear block P B and
    // shear-shear block S (shear order XY, YZ, XZ)
    const double A[3][3] = {{0.0, szz, syy}, {szz, 0.0, sxx}, {syy, sxx, 0.0}};
    const double B[3][3] = {{0.0, -2.0 * syz, 0.0}, {0.0, 0.0, -2.0 * sxz}, {-2.0 * sxy, 0.0, 0.0}};
    const double S[3][3] = {{-2.0 * szz, 2.0 * sxz, 2.0 * syz},
                            {2.0 * sxz, -2.0 * sxx, 2.0 * sxy},
                            {2.0 * syz, 2.0 * sxy, -2.0 * syy}};

    double PA[3][3], PAP[3][3], PB[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            PA[i][j] = 0.0;
            PB[i][j] = 0.0;
            for (int k = 0; k < 3; ++k)
            {
                PA[i][j] += P[i][k] * A[k][j];
                PB[i][j] += P[i][k] * B[k][j];
            }
        }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            PAP[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) PAP[i][j] += PA[i][k] * P[k][j];
        }

    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        for (int j = 0; j < VOIGTSIZE_3D; ++j)
        {
            double d2J2, d2J3;
            if (i < 3 && j < 3)
            {
                d2J2 = P[i][j];
                d2J3 = PAP[i][j];
            }
            else if (i < 3)
            {
                d2J2 = 0.0;
                d2J3 = PB[i][j - 3];
            }
            else if (j < 3)
            {
                d2J2 = 0.0;
                d2J3 = PB[j][i - 3];
            }
            else
            {
                d2J2 = (i == j) ? 2.0 : 0.0;
                d2J3 = S[i - 3][j - 3];
            }

            hessian[i * VOIGTSIZE_3D + j] = h_2 * d2J2 + h_3 * d2J3 + h_22 * dJ2[i] * dJ2[j] +
                                            h_23 * (dJ2[i] * dJ3[j] + dJ3[i] * dJ2[j]) +
                                            h_33 * dJ3[i] * dJ3[j];
        }
    }
}

/**
 * @brief cos(w - w_n) / (alpha C) for the deviatoric direction d(w) on the section, see
 * calculate_deviatoric_support_function.
 */
static double support_ratio(const double w, const double w_n, const MatsuokaNakaiConstants constants)
{
    double C, C1, C2;
    calculate_lode_function(-sin(3.0 * w), constants, &C, &C1, &C2);
    return cos(w - w_n) / (constants.alpha * C);
}

double calculate_deviatoric_support_function(const double j2, const double j3,
                                             const MatsuokaNakaiConstants constants)
{
    if (!(j2 > 0.0)) return 0.0;

    // The principal values 2 / sqrt(3) [sin(w + 2 pi / 3), sin(w), sin(w - 2 pi / 3)] describe the
    // deviatoric directions d(w) with J = 1 and sin(3 theta) = -sin(3 w); n = J_n d(w_n) up to the
    // order of its principal values. With s = d(w) / (alpha C) on the section,
    // s : n = 2 J_n cos(w - w_n) / (alpha C), which is maximised over w.
    const double J = sqrt(j2);
    double xi = 1.5 * sqrt(3.0) * j3 / (j2 * J);
    if (xi > 1.0) xi = 1.0;
    if (xi < -1.0) xi = -1.0;
    const double w_n = -asin(xi) / 3.0;

    // coarse search around w_n, then golden-section refinement around the best sample
    const int n_samples = 48;
    const double half_width = PI / 3.0;
    const double spacing = 2.0 * half_width / n_samples;
    double w_best = w_n, f_best = support_ratio(w_n, w_n, constants);
    for (int i = 0; i <= n_samples; ++i)
    {
        const double w = w_n - half_width + i * spacing;
        const double f = support_ratio(w, w_n, constants);
        if (f > f_best)
        {
            f_best = f;
            w_best = w;
        }
    }

    const double golden = 0.5 * (sqrt(5.0) - 1.0);
    double lo = w_best - spacing, hi = w_best + spacing;
    double w1 = hi - golden * (hi - lo), w2 = lo + golden * (hi - lo);
    double f1 = support_ratio(w1, w_n, constants), f2 = support_ratio(w2, w_n, constants);
    for (int iter = 0; iter < 40; ++iter)
    {
        if (f1 > f2)
        {
            hi = w2;
            w2 = w1;
            f2 = f1;
            w1 = hi - golden * (hi - lo);
            f1 = support_ratio(w1, w_n, constants);
        }
        else
        {
            lo = w1;
            w1 = w2;
            f1 = f2;
            w2 = lo + golden * (hi - lo);
            f2 = support_ratio(w2, w_n, constants);
        }
    }
    if (f1 > f_best) f_best = f1;
    if (f2 > f_best) f_best = f2;

    return 2.0 * J * f_best;
}

MatsuokaNakaiConstants calculate_matsuoka_nakai_constants(const double angle_rad, const double c)
{
    // closed form of the constants of Lagioia & Panteghini (2016) for k = (9 - sin^2) / (1 - sin^2),
    // A1 = (k - 3) / (k - 9), A2 = k / (k - 9): alpha = 2 / sqrt(3) sqrt(A1) M and
    // beta = A2 / A1^(3/2). Written in sin(angle), they are continuous at angle = 0, where the
    // surface becomes the von Mises cylinder J = 2 c / sqrt(3).
    const double s = sin(angle_rad);

    MatsuokaNakaiConstants constants;
    constants.M = 2.0 * sqrt(3.0) * s / (3.0 - s);
    constants.K = 2.0 * sqrt(3.0) * c * cos(angle_rad) / (3.0 - s);
    constants.alpha = 2.0 * sqrt(3.0 + s * s) / (3.0 - s);
    constants.beta = s * (9.0 - s * s) / pow(3.0 + s * s, 1.5);
    constants.gamma = 0.0;

    return constants;
}

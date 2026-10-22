#include <math.h>
#include <stddef.h>

#include "../globals.h"
#include "../strain_utils.h"
#include "../utils.h"
#include "../yield_surfaces/mohr_coulomb_surface.h"
#include "mohr_coulomb_return_mapping.h"

// tolerances relative to the magnitude of the trial stress
#define MC_YIELD_TOL 1.0e-12
#define MC_ORDER_TOL 1.0e-10
#define MC_LAMBDA_TOL 1.0e-10

/*
 * Return to the intersection of n (1 or 2) planes: the yield functions are linear in the stress,
 * so f_k(s_trial) = sum_l (a_k . D b_l) dLambda_l is solved directly.
 */
static int mc_return_to_planes(const double s_trial[3], const double D[9], int n,
                               const int planes[2][2], double sin_phi, double c_cos_phi,
                               double sin_psi, double s[3], double dlambda[2])
{
    double a[2][3], D_b[2][3], f[2], M[4];

    for (int k = 0; k < n; ++k)
    {
        double b[3];
        f[k] = mohr_coulomb_principal_function(s_trial, planes[k][0], planes[k][1], sin_phi,
                                               c_cos_phi, a[k]);
        mohr_coulomb_principal_gradient(planes[k][0], planes[k][1], sin_psi, b);
        matrix_vector_multiply(D, b, 3, D_b[k]);
    }
    for (int k = 0; k < n; ++k)
        for (int l = 0; l < n; ++l) M[2 * k + l] = vector_dot_product(a[k], D_b[l], 3);

    if (n == 1)
    {
        if (M[0] < SMALL_VALUE) return 0;
        dlambda[0] = f[0] / M[0];
    }
    else
    {
        double det = M[0] * M[3] - M[1] * M[2];
        if (fabs(det) < SMALL_VALUE * (M[0] * M[3] + fabs(M[1] * M[2]))) return 0;
        dlambda[0] = (M[3] * f[0] - M[1] * f[1]) / det;
        dlambda[1] = (M[0] * f[1] - M[2] * f[0]) / det;
    }

    for (int r = 0; r < 3; ++r)
    {
        s[r] = s_trial[r];
        for (int k = 0; k < n; ++k) s[r] -= dlambda[k] * D_b[k][r];
    }
    return 1;
}

static int mc_is_ordered(const double s[3], double tol)
{
    return (s[0] >= s[1] - tol && s[1] >= s[2] - tol) ? 1 : 0;
}

static void mc_set_return(MohrCoulombReturn* ret, MohrCoulombReturnType type, int n,
                          const int planes[2][2], const double dlambda[2])
{
    ret->type = type;
    ret->n_active = n;
    for (int k = 0; k < n; ++k)
    {
        ret->planes[k][0] = planes[k][0];
        ret->planes[k][1] = planes[k][1];
        ret->dlambda[k] = dlambda[k];
    }
}

int mohr_coulomb_return_mapping(const double s_trial[3], const double D[9], double sin_phi,
                                double cos_phi, double sin_psi, double c, double s[3],
                                MohrCoulombReturn* ret)
{
    static const int main_plane[2][2] = {{0, 2}, {0, 0}};
    static const int compression_edge[2][2] = {{0, 2}, {0, 1}};  // s2 = s3
    static const int extension_edge[2][2] = {{0, 2}, {1, 2}};    // s1 = s2

    double c_cos_phi = c * cos_phi;
    double scale = fmax(fabs(s_trial[0]), fabs(s_trial[2])) + c;
    double dlambda[2];

    ret->type = MC_RETURN_ELASTIC;
    ret->n_active = 0;
    copy_array(s_trial, 3, s);

    if (mohr_coulomb_principal_function(s_trial, 0, 2, sin_phi, c_cos_phi, NULL) <=
        MC_YIELD_TOL * scale)
        return 1;

    // main plane
    double s_plane[3];
    int plane_ok = mc_return_to_planes(s_trial, D, 1, main_plane, sin_phi, c_cos_phi, sin_psi,
                                       s_plane, dlambda);
    if (plane_ok && mc_is_ordered(s_plane, MC_ORDER_TOL * scale))
    {
        copy_array(s_plane, 3, s);
        mc_set_return(ret, MC_RETURN_PLANE, 1, main_plane, dlambda);
        return 1;
    }

    // edges, starting with the one of the principal stress pair whose order is lost most
    int compression_first = plane_ok ? (s_plane[2] - s_plane[1] >= s_plane[1] - s_plane[0]) : 1;
    for (int attempt = 0; attempt < 2; ++attempt)
    {
        int compression = (attempt == 0) ? compression_first : !compression_first;
        const int(*planes)[2] = compression ? compression_edge : extension_edge;
        double s_edge[3];
        if (!mc_return_to_planes(s_trial, D, 2, planes, sin_phi, c_cos_phi, sin_psi, s_edge,
                                 dlambda))
            continue;

        double lambda_tol = -MC_LAMBDA_TOL * (fabs(dlambda[0]) + fabs(dlambda[1]));
        if (dlambda[0] >= lambda_tol && dlambda[1] >= lambda_tol &&
            mc_is_ordered(s_edge, MC_ORDER_TOL * scale))
        {
            copy_array(s_edge, 3, s);
            mc_set_return(ret, compression ? MC_RETURN_EDGE_S2_EQ_S3 : MC_RETURN_EDGE_S1_EQ_S2, 2,
                          planes, dlambda);
            return 1;
        }
    }

    // apex
    if (sin_phi < SMALL_VALUE) return 0;
    double apex = -c * cos_phi / sin_phi;
    s[0] = s[1] = s[2] = apex;
    ret->type = MC_RETURN_APEX;
    ret->n_active = 0;
    return 1;
}

void mohr_coulomb_project_tangent(const MohrCoulombReturn* ret, double Q[3][3],
                                  const double De[VOIGTSIZE_3D * VOIGTSIZE_3D], double sin_phi,
                                  double sin_psi, double tangent[VOIGTSIZE_3D * VOIGTSIZE_3D])
{
    if (ret->type == MC_RETURN_ELASTIC) return;
    if (ret->type == MC_RETURN_APEX)
    {
        for (int i = 0; i < VOIGTSIZE_3D * VOIGTSIZE_3D; ++i) tangent[i] = 0.0;
        return;
    }

    int n = ret->n_active;
    double a6[2][VOIGTSIZE_3D], De_b[2][VOIGTSIZE_3D], AT_T[2][VOIGTSIZE_3D];
    double M[4], M_inv[4];

    for (int k = 0; k < n; ++k)
    {
        double a[3], b[3], b6[VOIGTSIZE_3D];
        mohr_coulomb_principal_gradient(ret->planes[k][0], ret->planes[k][1], sin_phi, a);
        mohr_coulomb_principal_gradient(ret->planes[k][0], ret->planes[k][1], sin_psi, b);
        calculate_strain_from_principal_system(a, Q, a6[k]);
        calculate_strain_from_principal_system(b, Q, b6);
        matrix_vector_multiply(De, b6, VOIGTSIZE_3D, De_b[k]);
    }

    for (int k = 0; k < n; ++k)
        for (int l = 0; l < n; ++l) M[2 * k + l] = vector_dot_product(a6[k], De_b[l], VOIGTSIZE_3D);
    if (n == 1)
    {
        M_inv[0] = 1.0 / M[0];
    }
    else
    {
        double det = M[0] * M[3] - M[1] * M[2];
        M_inv[0] = M[3] / det;
        M_inv[1] = -M[1] / det;
        M_inv[2] = -M[2] / det;
        M_inv[3] = M[0] / det;
    }

    // A^T T
    for (int l = 0; l < n; ++l)
        for (int j = 0; j < VOIGTSIZE_3D; ++j)
        {
            AT_T[l][j] = 0.0;
            for (int i = 0; i < VOIGTSIZE_3D; ++i)
                AT_T[l][j] += a6[l][i] * tangent[i * VOIGTSIZE_3D + j];
        }

    // T <- T - De N M^-1 A^T T
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
        for (int j = 0; j < VOIGTSIZE_3D; ++j)
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    tangent[i * VOIGTSIZE_3D + j] -= De_b[k][i] * M_inv[2 * k + l] * AT_T[l][j];
}

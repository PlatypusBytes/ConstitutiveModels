
#ifndef MATSUOKA_NAKAI_SURFACE_H
#define MATSUOKA_NAKAI_SURFACE_H



#ifdef __cplusplus
extern "C" {
#endif


void calculate_yield_function(const double p, const double theta, const double q,
                              const double c, const double phi_rad, const double alpha, const double beta,
                              const double gamma, const double K, const double M, double* f);

void calculate_yield_gradient(const double theta,
                               const double J,
                               const double* constants,
                               const double* dp_dsig, const double* dJ_dsig, const double* dtheta_dsig, double* grad);

void calculate_matsuoka_nakai_constants(const double angle_rad, const double c, double *alpha, double *beta, double *gamma, double *K, double *M);


#ifdef __cplusplus
}
#endif

#endif // MATSUOKA_NAKAI_SURFACE_H
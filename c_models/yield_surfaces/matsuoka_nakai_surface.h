
#ifndef MATSUOKA_NAKAI_SURFACE_H
#define MATSUOKA_NAKAI_SURFACE_H



#ifdef __cplusplus
extern "C" {
#endif


void calculate_yield_function(double p, double theta,double q,
                              double c, double phi_rad, double alpha, double beta, double gamma, double K, double M, double* f);

void calculate_yield_gradient(double theta,
                               double J,
                               double* constants,
                               double* dp_dsig, double* dJ_dsig, double* dtheta_dsig, double* grad);

void calculate_matsuoka_nakai_constants(double angle_rad, double c, double *alpha, double *beta, double *gamma, double *K, double *M);



#ifdef __cplusplus
}
#endif

#endif // MATSUOKA_NAKAI_SURFACE_H
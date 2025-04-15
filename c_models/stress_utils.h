#ifndef STRESS_UTILS_H
#define STRESS_UTILS_H

#ifdef __cplusplus
extern "C"
{
#endif

    void calculate_stress_invariants_3d(const double* stress, double* p, double* J, double* theta,
                                        double* j2, double* j3, double* s_dev);  // s_dev is output

    void calculate_stress_invariants_derivatives_3d(const double J, const double* s_dev,
                                                    const double j2, const double j3,
                                                    double* dp_dsig, double* dJ_dsig,
                                                    double* dtheta_dsig);

#ifdef __cplusplus
}
#endif

#endif  // STRESS_UTILS_H

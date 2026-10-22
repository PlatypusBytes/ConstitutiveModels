#include "rowe_dilatancy.h"

double rowe_mobilised_sin_psi(double sin_phi_m, double sin_phi_cv)
{
    return (sin_phi_m - sin_phi_cv) / (1.0 - sin_phi_m * sin_phi_cv);
}

double rowe_critical_state_sin_phi(double sin_phi, double sin_psi)
{
    return (sin_phi - sin_psi) / (1.0 - sin_phi * sin_psi);
}

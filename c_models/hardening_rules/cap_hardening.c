#include "cap_hardening.h"

double cap_hardening_modulus(double Ks, double K_ratio)
{
    return Ks / (K_ratio - 1.0);
}

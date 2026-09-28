"""
Python counterpart of c_models/hardening_rules/hyperbolic_shear_hardening.c and cap_hardening.c.
"""


class HardeningLaws:

    @staticmethod
    def hyperbolic_plastic_shear_strain(q, qa, Ei, Eur):
        """Plastic shear strain gamma_p on the hyperbola at deviator q (Eq. 8 with f = 0)."""
        return 2.0 / Ei * q / (1.0 - q / qa) - 2.0 * q / Eur

    @staticmethod
    def hyperbolic_initial_stiffness(E50, Rf):
        """Initial stiffness Ei such that E50 is the secant stiffness at q = qf / 2 (Sec. 2.1)."""
        return 2.0 * E50 / (2.0 - Rf)

    @staticmethod
    def cap_hardening_modulus(Ks, K_ratio):
        """Cap hardening modulus H = Ks Kc / (Ks - Kc) = Ks / (Ks/Kc - 1) (Eq. 32)."""
        return Ks / (K_ratio - 1.0)

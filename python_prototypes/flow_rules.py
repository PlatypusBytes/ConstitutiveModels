"""
Python counterpart of c_models/flow_rules/rowe_dilatancy.c.
"""


class FlowRules:

    @staticmethod
    def rowe_mobilised_sin_psi(sin_phi_m, sin_phi_cv):
        """Mobilised dilatancy of Rowe's stress-dilatancy theory (Eq. 11)."""
        return (sin_phi_m - sin_phi_cv) / (1.0 - sin_phi_m * sin_phi_cv)

    @staticmethod
    def rowe_critical_state_sin_phi(sin_phi, sin_psi):
        """Critical state friction angle from the friction and dilation angle at failure (Eq. 13)."""
        return (sin_phi - sin_psi) / (1.0 - sin_phi * sin_psi)


from python_models.solving_algorithms.brute_force import BruteForceAlgorithm
from python_models.utils import StressUtils

import numpy as np


class YieldSurfaceCorrectionConstitutiveModel:
    """
    For general elastoplastic models

    sloan et al 2001 page 135
    """

    def __init__(self, yield_function, flow_rule, elasticity_model, solving_algorithm=BruteForceAlgorithm(), hardening_rule= None, plastic_multiplier_rule= None):

        self.yield_function = yield_function
        self.flow_rule = flow_rule
        self.elasticity_model = elasticity_model
        self.solving_algorithm = solving_algorithm

        self.stress_utils = StressUtils()

        self.solving_algorithm = None
        self.hardening_rule = None
        self.plastic_multiplier_rule = None

        self.D_el = None
        self.D_tot = None


    def calculate_plastic_potential_multiplier(self, strain_vector):
        """
        Calculate the increment of the plastic multiplier

        :return:
        """

        self.stress_utils.update_derivatives()

        self.yield_function.calculate_derivative()
        derivative_yield_function = self.yield_function.derivative

        self.flow_rule.calculate_derivative()
        derivative_flow_rule = self.flow_rule.derivative

        # standard approach
        # delta_lambda = np.inner(np.inner(derivative_yield_function,self.D_el),strain_vector) / (
        #     np.inner(np.inner(derivative_yield_function,self.D_el),derivative_flow_rule) )

        # euler backward approach
        delta_lambda = self.yield_function.yield_residual / np.inner(np.inner(derivative_yield_function,self.D_el),derivative_flow_rule)

        if np.isnan(delta_lambda):
            raise ValueError('delta lambda is nan')

        return delta_lambda



    def update_stress(self,strain_vector):

        delta_lambda = self.calculate_plastic_potential_multiplier(strain_vector)

        # standard approach, trial stress - delta_lambda * D_el.dot(flow_rule.derivative)
        stress = self.stress_utils.stress_vector - delta_lambda *self.D_el.dot(self.flow_rule.derivative)

        return stress

    def consistency_condition(self, stress):

        # stress = self.stress_utils.stress_vector - factor * self.D_el.dot(self.flow_rule.derivative)
        self.stress_utils.update_stress(stress)
        self.stress_utils.update_invariants()

        self.yield_function.calculate_yield_function()

        return self.yield_function.yield_residual

    def solve_1(self, strain_vector):

        convergence_tol = 1e-6
        max_internal_iteration = 100
        i = 0
        y0 = self.yield_function.yield_residual
        # if the yield function is greater than zero, then plasticity is activated
        while y0 > convergence_tol and i < max_internal_iteration:
            delta_lambda = self.calculate_plastic_potential_multiplier(strain_vector)
            stress = self.stress_utils.stress_vector - delta_lambda * self.D_el.dot(self.flow_rule.derivative)

            y0 = self.consistency_condition(stress)

            i += 1

        if i == max_internal_iteration:
            raise ValueError('No convergence max internal iteration reached')

    def solve_pegasus(self, strain_vector):
        convergence_tol = 1e-6
        max_internal_iteration = 100
        i = 0
        F_new = self.yield_function.yield_residual
        alpha_0 = 0
        alpha_1 = 1

        delta_sig_trial = self.D_el.dot(strain_vector)
        F0 = self.consistency_condition(self.original_stress + alpha_0 * delta_sig_trial)
        F1 = self.consistency_condition(self.original_stress + alpha_1 * delta_sig_trial)
        alpha_new =1
        # if the yield function is greater than zero, then plasticity is activated
        while F_new > convergence_tol and i < max_internal_iteration:
            alpha_new = alpha_1 - F1 * (alpha_1 - alpha_0) / (F1 - F0)
            delta_sig_trial = self.D_el.dot(strain_vector)
            F_new = self.consistency_condition(self.original_stress + alpha_new * delta_sig_trial)
            # if F_new is opposite sign to F0
            if F_new * F0 < 0:
                alpha_1 = alpha_new
                F1 = F_new
            else:
                F1 = F1*F0/(F0*F_new)
                alpha_0 = alpha_new
                F0 = F_new

            i+=1

        a=1+1
        # self.consistency_condition(self.original_stress + alpha_new * delta_sig_trial)



    def main(self, stress_vector, strain_vector):
        """

        :param stress_vector:
        :param strain_vector:
        :param material_parameters:
        :return:
        """

        self.original_stress = np.copy(stress_vector)

        self.stress_utils.update_stress(stress_vector)

        self.D_el = self.elasticity_model.calculate_elastic_matrix()

        trial_delta_stress = self.D_el.dot(strain_vector)
        trial_stress = self.stress_utils.stress_vector + trial_delta_stress

        self.stress_utils.update_stress(trial_stress)
        self.stress_utils.update_invariants()

        self.yield_function.stress_utils = self.stress_utils
        self.yield_function.initialize()
        self.yield_function.calculate_yield_function()

        convergence_tol = 1e-6
        max_internal_iteration = 100
        i = 0

        if self.yield_function.yield_residual > 0:
            self.flow_rule.stress_utils = self.stress_utils
            self.flow_rule.initialize()

        # self.solve_1(strain_vector)
        self.solve_pegasus(strain_vector)

        # y0= self.yield_function.yield_residual
        # # if the yield function is greater than zero, then plasticity is activated
        # while y0 > convergence_tol and i < max_internal_iteration:
        #
        #     delta_lambda = self.calculate_plastic_potential_multiplier(strain_vector)
        #     stress = self.stress_utils.stress_vector - delta_lambda * self.D_el.dot(self.flow_rule.derivative)
        #     y0 = self.consistency_condition(stress)
        #
        #     # # elastoplastic stress update
        #     # stress = self.update_stress(strain_vector)
        #     # self.stress_utils.update_stress(stress)
        #     # self.stress_utils.update_invariants()
        #     #
        #     # self.yield_function.calculate_yield_function()
        #     i += 1

        if i == max_internal_iteration:
            raise ValueError('No convergence max internal iteration reached')


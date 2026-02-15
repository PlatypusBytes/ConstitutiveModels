
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

        self.A = 0 # hardening parameter -dfdh dh/dlambda
        self.B = 0 # hardening parameter -A/ (df/dh)
        self.H = 0 # initial hardening parameter


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
        delta_lambda = self.yield_function.yield_residual /(self.A + np.inner(np.inner(derivative_yield_function,self.D_el),derivative_flow_rule))

        if np.isnan(delta_lambda):
            raise ValueError('delta lambda is nan')

        return delta_lambda



    # def update_stress(self,strain_vector):
    #
    #     delta_lambda = self.calculate_plastic_potential_multiplier(strain_vector)
    #
    #     # standard approach, trial stress - delta_lambda * D_el.dot(flow_rule.derivative)
    #     stress = self.stress_utils.stress_vector - delta_lambda *self.D_el.dot(self.flow_rule.derivative)
    #
    #     return stress


    def solve(self, strain_vector):

        convergence_tol = 1e-6
        max_internal_iteration = 100
        i = 0
        y0 = self.yield_function.yield_residual
        H0 = self.H
        # if the yield function is greater than zero, then plasticity is activated
        uncorrected_stress = np.copy(self.stress_utils.stress_vector)
        while abs(y0) > convergence_tol and i < max_internal_iteration:
            delta_lambda = self.calculate_plastic_potential_multiplier(strain_vector)
            corrected_stress = uncorrected_stress - delta_lambda * self.D_el.dot(self.flow_rule.derivative)
            self.H = self.H + self.B * delta_lambda

            self.update_stress(corrected_stress)

            self.yield_function.calculate_yield_function()
            new_y0 = self.yield_function.yield_residual

            if abs(new_y0) > abs(y0):
                self.update_stress(uncorrected_stress)
                self.stress_utils.update_derivatives()
                self.yield_function.calculate_derivative()
                delta_lambda = y0 / np.inner(self.yield_function.derivative, self.yield_function.derivative)
                corrected_stress = uncorrected_stress - delta_lambda * self.yield_function.derivative
                self.update_stress(corrected_stress)
                self.H = H0

                self.yield_function.calculate_yield_function()
                new_y0 = self.yield_function.yield_residual

            uncorrected_stress = np.copy(corrected_stress)
            H0 = self.H
            y0 = new_y0

            i+=1

        if i == max_internal_iteration:
            raise ValueError('No convergence max internal iteration reached')


    def update_stress(self, stress_vector):
        self.stress_utils.update_stress(stress_vector)
        self.stress_utils.update_invariants()

    def main(self, stress_vector, strain_vector):
        """

        :param stress_vector:
        :param strain_vector:
        :param material_parameters:
        :return:
        """
        self.yield_function.initialize()
        self.yield_function.stress_utils = self.stress_utils


        self.original_stress = np.copy(stress_vector)
        # self.stress_utils.update_stress(stress_vector)
        self.D_el = self.elasticity_model.calculate_elastic_matrix()

        trial_delta_stress = self.D_el.dot(strain_vector)
        trial_stress = self.original_stress + trial_delta_stress

        self.update_stress(trial_stress)

        self.yield_function.calculate_yield_function()

        # self.stress_utils.update_stress(trial_stress)
        # self.stress_utils.update_invariants()
        #
        # self.yield_function.stress_utils = self.stress_utils
        # # self.yield_function.initialize()
        # self.yield_function.calculate_yield_function()

        convergence_tol = 1e-6
        max_internal_iteration = 100
        i = 0

        if self.yield_function.yield_residual > 0:
            self.flow_rule.stress_utils = self.stress_utils
            self.flow_rule.initialize()

            self.solve(strain_vector)
        # self.solve_pegasus(strain_vector)

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


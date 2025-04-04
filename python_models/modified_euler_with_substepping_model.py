
from python_models.solving_algorithms.brute_force import BruteForceAlgorithm
from python_models.utils import StressUtils

import numpy as np


class ModifiedEulerWithSubsteppingModel:
    """
    For critical state models

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
        stress = np.copy(stress_vector)
        # self.stress_utils.update_stress(stress_vector)
        self.D_el = self.elasticity_model.calculate_elastic_matrix()

        trial_delta_stress = self.D_el.dot(strain_vector)
        trial_stress = self.original_stress + trial_delta_stress

        self.update_stress(trial_stress)
        self.yield_function.calculate_yield_function()
        y_el = self.yield_function.yield_residual

        tol =1e-6
        alpha = 0
        if self.yield_function.yield_residual > 0:



            self.flow_rule.stress_utils = self.stress_utils
            self.flow_rule.initialize()

            # check if there is a transition from elastic to plastic
            self.update_stress(self.original_stress)
            self.yield_function.calculate_yield_function()
            y0 = self.yield_function.yield_residual

            if y0 < -tol and y_el > tol:
                #pegasus algorithm
                alpha = 0.5
                pass

            if abs(y0) < tol and y_el > tol:
                # check elastoplastic unloading
                self.solve(strain_vector)


            # update stress
            self.elasticity_model.calculate_elastic_matrix()
            stress = stress + alpha * self.elasticity_model.elastic_matrix.dot(strain_vector)
            strain_vector = (1-alpha) * strain_vector
            T= 0
            dT = 1
            while T<1:
                for i in range(2):
                    delta_stress_el_i = self.elasticity_model.elastic_matrix.dot(strain_vector)
                    delta_stress_i = delta_stress_el_i

                    self.update_stress(stress + delta_stress_i)

                    self.yield_function.calculate_derivative()
                    dfdsigma = self.yield_function.derivative
                    self.flow_rule.calculate_derivative()
                    dgdsigma = self.flow_rule.derivative

                    Bi = 0
                    dfdh = 0
                    Ai = -dfdh * Bi


                    delta_lambda = max(0, np.inner(dfdsigma, delta_stress_i) /(Ai+ np.inner(dfdsigma, dgdsigma)))
                    delta_H = delta_lambda * Bi


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



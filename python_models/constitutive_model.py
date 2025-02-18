
from models.solving_algorithms.brute_force import BruteForceAlgorithm
from models.utils import StressUtils

import numpy as np

class BaseConstitutiveModel:

    def __init__(self, yield_function, flow_rule, elasticity_equation, solving_algorithm=BruteForceAlgorithm(), hardening_rule= None, plastic_multiplier_rule= None):

        self.yield_function = yield_function
        self.flow_rule = flow_rule
        self.elasticity_equation = elasticity_equation
        self.solving_algorithm = solving_algorithm

        self.stress_utils = StressUtils()

        self.solving_algorithm = None
        self.hardening_rule = None
        self.plastic_multiplier_rule = None

        self.D_el = None
        self.D_tot = None


    def calculate_plastic_potential_multiplier(self):
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
        # delta_lambda1 = np.inner(np.inner(derivative_yield_function,D),self.strain_vector) / (
        #     np.inner(np.inner(derivative_yield_function,D),derivative_dilatancy_function) )

        # euler backward approach
        delta_lambda = self.yield_function.yield_residual / np.inner(np.inner(derivative_yield_function,self.D_el),derivative_flow_rule)

        if np.isnan(delta_lambda):
            raise ValueError('delta lambda is nan')

        return delta_lambda



    def update_stress(self):

        delta_lambda = self.calculate_plastic_potential_multiplier()

        # standard approach, trial stress - delta_lambda * D_el.dot(flow_rule.derivative)
        stress = self.stress_utils.stress_vector - delta_lambda *self.D_el.dot(self.flow_rule.derivative)


        return stress


    def main(self, stress_vector, strain_vector):
        """

        :param stress_vector:
        :param strain_vector:
        :param material_parameters:
        :return:
        """

        self.stress_utils.update_stress(stress_vector)

        self.D_el = self.elasticity_equation.calculate_elastic_matrix()

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


        # if the yield function is greater than zero, then plasticity is activated
        while self.yield_function.yield_residual > convergence_tol and i < max_internal_iteration:

            # elastoplastic stress update
            stress = self.update_stress()
            self.stress_utils.update_stress(stress)
            self.stress_utils.update_invariants()

            self.yield_function.calculate_yield_function()
            i += 1

        if i == max_internal_iteration:
            raise ValueError('No convergence max internal iteration reached')



if __name__ == '__main__':

    from models.incr_driver import IncrDriver
    from models.yield_surfaces.matsuoka_nakai import MatsuokaNakai
    from models.elasticity_models.hookes_law import HooksLaw

    yield_function = MatsuokaNakai({'angle': 30, 'cohesion': 0.00})
    flow_rule = MatsuokaNakai({'angle': 0, 'cohesion': 0.00})

    G = 2000
    nu = 0.33
    E = 2 * G * (1 + nu)
    elasticity_model = HooksLaw({'young_modulus': E, 'poison_ratio': nu})
    elasticity_model.calculate_elastic_matrix()
    #
    # params = {'young_modulus': E, 'poison_ratio': nu, 'phi': 30, 'psi': 1e-5, 'c': 0.00, 'tensile_cutoff': 10000}
    #
    stress_vector = np.array([-1, -1, -1, 0, 0, 0])
    strain_increment = -0.0001
    # d_stran = -0.0001

    strain = np.zeros(6)
    d_strain = np.zeros(6)
    #
    # delta_strain = np.array([0,0,-0.0001,0,0,0])
    incr_stress_input = np.array([0, 0, 0, 0, 0, 0])
    control_type = [1, 1, 0, 1, 1, 1]

    orig_stress_vector = np.copy(stress_vector)
    all_stress_vectors = [stress_vector]
    all_strain_vectors = [np.zeros(6)]

    delta_stress = np.zeros(6)
    strain_vector = np.zeros(6)


    const_model = BaseConstitutiveModel(yield_function, flow_rule, elasticity_model)
    incr_driver = IncrDriver()

    for t in range(10):

        delta_strain = np.array([0, 0, strain_increment, 0, 0, 0])
        correction_delta_strain = np.zeros(6)

        approx_delta_stress = np.zeros(6)
        old_d_stress_vector = np.copy(stress_vector)

        max_iter = 100
        for i in range(max_iter + 1):

            u_d_stress = np.zeros(6) - approx_delta_stress
            incr_driver.calculate_time_step(correction_delta_strain, u_d_stress, control_type, elasticity_model.elastic_matrix)
            delta_strain = delta_strain + correction_delta_strain

            const_model.main(stress_vector, delta_strain)

            if i < max_iter:
                approx_delta_stress = const_model.stress_utils.stress_vector - old_d_stress_vector
                stress_vector = np.copy(old_d_stress_vector)
            else:
                strain_vector = strain_vector + delta_strain
                stress_vector = const_model.stress_utils.stress_vector

        all_stress_vectors.append(stress_vector)
        all_strain_vectors.append(strain_vector)

    a = 1 + 1
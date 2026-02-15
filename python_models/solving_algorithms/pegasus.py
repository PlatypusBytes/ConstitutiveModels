import numpy as np

class Pegasus(SolvingAlgorithm):
    def __init__(self, consistency_equation, elasticity_model):
        super().__init__()
        self.consistency_equation = consistency_equation
        self.elasticity_model = elasticity_model
        self.alpha_0 = 0
        self.alpha_1 = 1


    def solve(self, strain_vector):

        # calculate the elastic matrix 0 and 1, for now its constant
        D_el_0 = np.copy(self.elasticity_model.elastic_matrix)
        D_el_1 = np.copy(self.elasticity_model.elastic_matrix)

        delta_sig_0 = self.alpha_0 * D_el_0.dot(strain_vector)
        delta_sig_1 = self.alpha_1 * D_el_1.dot(strain_vector)

        F0 = self.consistency_equation(strain_vector)
        F1 = self.consistency_equation(strain_vector)

        max_its = 100
        tol = 1e-6
        i = 0
        while self.consistency_equation.yield_residual > tol and i < max_its:

            self.consistency_equation.update_stress()
            self.consistency_equation.update_invariants()
            self.consistency_equation.calculate_yield_function()
            i += 1

        return self.solver(instance, self.max_iter, self.max_time, self.seed)

from python_models.solving_algorithms.solving_algorithm_abc import SolvingAlgorithmABC

class BruteForceAlgorithm(SolvingAlgorithmABC):
    def __init__(self ):
        super().__init__()

    #
    # def solve(self):
    #
    #     # if the yield function is greater than zero, then plasticity is activated
    #     while self.yield_function.yield_residual > convergence_tol and i < max_internal_iteration:
    #         # elastoplastic stress update
    #         stress = self.update_stress(strain_vector)
    #         self.stress_utils.update_stress(stress)
    #         self.stress_utils.update_invariants()
    #
    #         self.yield_function.calculate_yield_function()
    #         i += 1
    #
    #     if i == max_internal_iteration:
    #         raise ValueError('No convergence max internal iteration reached')
    #
    #     # Solve the problem using brute force
    #     pass

from models.elasticity_models.base_law import BaseLaw

import numpy as np

class HooksLaw(BaseLaw):
    def __init__(self, params):
        super().__init__()
        self.params = params

        self.elastic_matrix = None



    def calculate_elastic_matrix(self) -> np.ndarray:
        """
        Calculate the elastic matrix for the given material parameters

        :param params: dict: Material parameters
        :return D: Elastic matrix
        """

        E = self.params['young_modulus']
        nu = self.params['poison_ratio']

        # Elastic matrix
        D = np.zeros((6,6))

        D[0,0] = D[1,1] = D[2,2] = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
        D[3,3] = D[4,4] = D[5,5] = E / (2 * (1 + nu))
        D[0,1] = D[0,2] = D[1,0] = D[1,2] = D[2,0] = D[2,1] = E * nu / ((1 + nu) * (1 - 2 * nu))

        self.elastic_matrix = D

        return D


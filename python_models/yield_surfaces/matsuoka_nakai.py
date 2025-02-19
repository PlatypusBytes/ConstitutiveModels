from python_models.yield_surfaces.general_yield_surface import GeneralYieldSurface

import numpy as np



class MatsuokaNakai(GeneralYieldSurface):

    def __init__(self,params):
        super().__init__()

        self.angle = params['angle'] # in degrees
        self.cohesion = params['cohesion']


    def initialize(self):
        self.calculate_constants()

    def calculate_constants(self):
        # if angle is zero, return tresca constants
        if np.isclose(self.angle, 0):
            self.M = 0
            self.K = 0
            self.alpha = 1 / (np.cos(np.pi / 6))
            self.beta = 0.9999
            self.gamma = 1

        # else Matsuaoka Nakai constants
        else:

            self.M = 1 / np.sqrt(3) * 6 * np.sin(np.radians(self.angle)) / (3 - np.sin(np.radians(self.angle)))
            self.K = self.cohesion / np.tan(np.radians(self.angle))
            k_mn = (9 - np.sin(np.radians(self.angle)) ** 2) / (1 - np.sin(np.radians(self.angle)) ** 2)
            A1 = (k_mn - 3) / (k_mn - 9)
            A2 = k_mn / (k_mn - 9)

            self.alpha = 2 / np.sqrt(3) * np.sqrt(A1) * self.M
            self.beta = A2 / (A1 ** (3 / 2))
            self.gamma = 0

    def calculate_yield_function(self):
        super().calculate_yield_function()

    def calculate_derivative(self):
        super().calculate_derivative()


from abc import abstractmethod

import numpy as np


class GeneralYieldSurface():
    def __init__(self):
        self.K = None
        self.M = None
        self.alpha = None
        self.beta = None
        self.gamma = None

        self.yield_residual = None
        self.derivative = None

        self.stress_utils = None
        pass

    @abstractmethod
    def initialize(self):
        # raise error because this is an abstract method
        raise NotImplemented('abstract method of initialize is called')

    @abstractmethod
    def calculate_constants(self, *args):
        # raise error because this is an abstract method
        raise NotImplemented('abstract method of calculate_constants is called')

    def calculate_yield_function(self):

        J = np.sqrt(self.stress_utils.second_deviatoric_invariant)

        self.yield_residual = (-self.K + self.M * self.stress_utils.mean_stress) + J * self.alpha * np.cos(
            np.acos(self.beta * np.sin(3 * self.stress_utils.lode_angle)) / 3 - self.gamma * np.pi / 6)


    def calculate_derivative(self):
        J2 = self.stress_utils.second_deviatoric_invariant
        J = np.sqrt(J2)
        theta = self.stress_utils.lode_angle

        dmean_dsigma = self.stress_utils.derivative_mean_stress
        dJ_dsigma = self.stress_utils.derivative_second_deviatoric_invariant / (2 * np.sqrt(J2))
        dtheta_dsigma = self.stress_utils.derivative_lode_angle

        # angle_term = np.acos(self.beta_phi*np.sin(3*theta))/3
        angle_term = np.acos(self.beta * np.sin(3 * theta)) / 3 - self.gamma * np.pi / 6

        dfdsigma = (self.M * dmean_dsigma + self.alpha *
                    (dJ_dsigma * (np.cos(angle_term)) + J * np.sin(angle_term) * self.beta * np.cos(
                        3 * theta) * dtheta_dsigma) / np.sqrt(1 - self.beta ** 2 * np.sin(3 * theta) ** 2))

        if np.isnan(dfdsigma).any():
            raise ValueError('dfdsigma is nan')

        self.derivative = dfdsigma

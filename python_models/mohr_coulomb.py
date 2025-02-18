
from abc import abstractmethod

import numpy as np

from models.utils import StressUtils
from models.incr_driver import IncrDriver

class GeneralYieldFunction():
    def __init__(self):
        self.K = None
        self.M = None
        self.alpha = None
        self.beta = None
        self.gamma = None

        self.yield_residual = None
        pass

    @abstractmethod
    def calculate_constants(self, *args):
        # raise error because this is an abstract method
        raise NotImplemented('abstract method of calculate_constants is called')

    def calculate_yield_function(self, mean_stress, second_deviatoric_invariant, lode_angle):

        J = np.sqrt(second_deviatoric_invariant)

        self.yield_residual = (-self.K + self.M * mean_stress) + J * self.alpha * np.cos(
            np.acos(self.beta * np.sin(3 * lode_angle)) / 3 - self.gamma * np.pi / 6)



class MatsuokaNakai(GeneralYieldFunction):
    def __init__(self, stress_tensor, strain_vector, params ,t):
        super().__init__()

        self.strain_vector = strain_vector
        self.stress_utils = StressUtils(stress_tensor)

        self.phi = params['phi']
        self.psi = params['psi']
        self.c = params['c']
        self.poison_ratio = params['poison_ratio']
        self.tensile_cutoff = params['tensile_cutoff']

        self.E = params['E']
        self.time_step = t

        self.D_el = None
        self.D_pl = None


    # def calculate_constants(self, angle, cohesion):
    #
    #     # if angle is zero, return tresca constants
    #     if np.isclose(angle,0):
    #         self.M = 0
    #         self.K = 0
    #         self.alpha = 1/(np.cos(np.pi/6))
    #         self.beta = 0.9999
    #         self.gamma = 1
    #         # return M, K, alpha, beta, gamma
    #
    #     # else matsuaoka nakai constants
    #     else:
    #
    #
    #         self.M = 1/np.sqrt(3) * 6* np.sin(np.radians(angle)) / (3-np.sin(np.radians(angle)))
    #         self.K = cohesion / np.tan(np.radians(angle))
    #         k_mn = (9-np.sin(np.radians(angle))**2)/(1-np.sin(np.radians(angle))**2)
    #         A1 = (k_mn -3)/(k_mn -9)
    #         A2 = k_mn/(k_mn -9)
    #
    #         self.alpha = 2/np.sqrt(3) * np.sqrt(A1) * self.M
    #         self.beta = A2/(A1**(3/2))
    #         self.gamma = 0
    #
    #
    # def calculate_friction_constants(self):
    #     self.M_phi, self.K_phi, self.alpha_phi, self.beta_phi, self.gamma_phi = (
    #         self.calculate_constants(self.phi, self.c))
    #
    # def calculate_dilatancy_constants(self):
    #     self.M_psi, self.K_psi, self.alpha_psi, self.beta_psi, self.gamma_psi = (
    #         self.calculate_constants(self.psi, self.c))

    # def calculate_yield_surface(self):
    #
    #     J2 = self.stress_utils.get_second_deviatoric_invariant()
    #     J = np.sqrt(J2)
    #     p = self.stress_utils.get_mean_stress()
    #     theta = self.stress_utils.get_lode_angle()
    #
    #     self.calculate_friction_constants()
    #
    #     yield_function = (-self.K_phi + self.M_phi*p) + J * self.alpha_phi * np.cos(np.acos(self.beta_phi*np.sin(3*theta))/3 -self.gamma_phi*np.pi/6)
    #
    #     return yield_function
    #
    # def calculate_tensile_yield_function(self):
    #     J2 = self.stress_utils.get_second_deviatoric_invariant()
    #     I1 = self.stress_utils.get_first_invariant()
    #
    #     theta = self.stress_utils.get_lode_angle()
    #
    #     yield_function = -(self.tensile_cutoff - (I1/3 + 2*np.sqrt(J2/3)*np.cos(theta-4*np.pi/3)))
    #
    #     return yield_function
    #
    # def calculate_derivative_yield_function(self):
    #     J2 = self.stress_utils.get_second_deviatoric_invariant()
    #     J = np.sqrt(J2)
    #     theta = self.stress_utils.get_lode_angle()
    #
    #     dmean_dsigma = self.stress_utils.get_derivative_mean_stress()
    #     dJ_dsigma = self.stress_utils.get_derivative_second_deviatoric_invariant()/(2*np.sqrt(J2))
    #     dtheta_dsigma = self.stress_utils.get_derivative_lode_angle()
    #
    #     # angle_term = np.acos(self.beta_phi*np.sin(3*theta))/3
    #     angle_term = np.acos(self.beta_phi * np.sin(3 * theta)) / 3 - self.gamma_phi * np.pi / 6
    #
    #     dfdsigma = (self.M_phi * dmean_dsigma + self.alpha_phi*
    #                 (dJ_dsigma * (np.cos(angle_term)) + J * np.sin(angle_term)*self.beta_phi*np.cos(3*theta) * dtheta_dsigma)/np.sqrt(1-self.beta_phi**2*np.sin(3*theta)**2))
    #
    #     if np.isnan(dfdsigma).any():
    #         raise ValueError('dfdsigma is nan')
    #
    #     return dfdsigma

    def calculate_dilatancy_function(self):
        J2 = self.stress_utils.get_second_deviatoric_invariant()
        J = np.sqrt(J2)
        p = self.stress_utils.get_mean_stress()
        theta = self.stress_utils.get_lode_angle()

        self.calculate_dilatancy_constants()

        dilatancy_function = (-self.K_psi + self.M_psi*p) + J * self.alpha_psi * np.cos(np.acos(self.beta_psi*np.sin(3*theta))/3 -self.gamma_psi*np.pi/6)
        return dilatancy_function

    def calculate_derivative_dilatancy_function(self):
        J2 = self.stress_utils.get_second_deviatoric_invariant()
        J = np.sqrt(J2)
        theta = self.stress_utils.get_lode_angle()

        dmean_dsigma = self.stress_utils.get_derivative_mean_stress()
        dJ_dsigma = self.stress_utils.get_derivative_second_deviatoric_invariant()/(2*np.sqrt(J2))
        dtheta_dsigma = self.stress_utils.get_derivative_lode_angle()

        # angle_term = np.acos(self.beta_psi*np.sin(3*theta))/3
        angle_term = np.acos(self.beta_psi * np.sin(3 * theta)) / 3 - self.gamma_psi * np.pi / 6

        dgdsigma = (self.M_psi * dmean_dsigma + self.alpha_psi*
                   (dJ_dsigma * (np.cos(angle_term)) + J * self.beta_psi* np.sin(angle_term)*np.cos(3*theta) * dtheta_dsigma)/
                   np.sqrt(1-self.beta_psi**2*np.sin(3*theta)**2))

        return dgdsigma

    # def calculate_elastic_matrix(self):
    #     D = np.zeros((6,6))
    #     young_modulus = self.E
    #     poisson_ratio = self.poison_ratio
    #     D[0,0] = D[1,1] = D[2,2] = young_modulus * (1 - poisson_ratio) / ((1 + poisson_ratio) * (1 - 2 * poisson_ratio))
    #     D[3,3] = D[4,4] = D[5,5] = young_modulus / (2 * (1 + poisson_ratio))
    #     D[0,1] = D[0,2] = D[1,0] = D[1,2] = D[2,0] = D[2,1] = young_modulus * poisson_ratio / ((1 + poisson_ratio) * (1 - 2 * poisson_ratio))
    #
    #     self.D_el = D


    def calculate_delta_lambda(self):
        """
        Calculate the increment of the plastic multiplier

        :return:
        """

        derivative_yield_function = self.calculate_derivative_yield_function()
        derivative_dilatancy_function = self.calculate_derivative_dilatancy_function()


        # standard approach
        # delta_lambda1 = np.inner(np.inner(derivative_yield_function,D),self.strain_vector) / (
        #     np.inner(np.inner(derivative_yield_function,D),derivative_dilatancy_function) )

        # euler backward approach
        delta_lambda = self.calculate_yield_function() / np.inner(np.inner(derivative_yield_function,self.D_el),derivative_dilatancy_function)

        if np.isnan(delta_lambda):
            raise ValueError('delta lambda is nan')

        return delta_lambda


    def calculate_elasto_plastic_matrix(self):

        # delta_lambda = self.calculate_delta_lambda()
        dgdsigma = self.calculate_derivative_dilatancy_function()
        dfdsigma = self.calculate_derivative_yield_function()
        D_ep = self.D_el - self.D_el.dot(dgdsigma)[:,None].dot(dfdsigma[None, :].dot(self.D_el))/(dfdsigma[None,:].dot(self.D_el).dot(dgdsigma))
        # D_bar = D - D.dot(dgdsigma).dot(dfdsigma[:,None].dot(D))/(dfdsigma[None,:].dot(D).dot(dgdsigma))
        return D_ep

    def update_stress(self):

        delta_lambda = self.calculate_delta_lambda()

        # standard approach
        stress = self.stress_utils.stress_vector - delta_lambda *self.D_el.dot(self.calculate_derivative_dilatancy_function())


        return stress

    def main(self):

        self.calculate_elastic_matrix()

        delta_sig = self.D_el.dot(self.strain_vector)
        self.stress_utils.stress_vector = self.stress_utils.stress_vector + delta_sig

        self.stress_utils.calculate_stress_tensor()

        # calculate yield function
        self.calculate_friction_constants()
        self.calculate_dilatancy_constants()
        f0 = self.calculate_yield_function()

        # f_tensile = self.calculate_tensile_yield_function()

        # if f_tensile > 0.0:
        #     a=1+1

        convergence_tol = 1e-6
        max_internal_iteration = 100
        i = 0
        while f0 > convergence_tol and i < max_internal_iteration:

            #elastoplastic
            stress = self.update_stress()

            self.stress_utils.stress_vector = stress
            self.stress_utils.calculate_stress_tensor()

            f0 = self.calculate_yield_function()
            i += 1

        if i == max_internal_iteration:
            raise ValueError('No convergence max internal iteration reached')

stress_tensor = np.array([[-1, 0, 0], [0,-1, 0], [0, 0, -1]])

G =2000
nu = 0.33
E = 2 * G * (1 + nu)

params = {'E':E,'poison_ratio':nu,'phi': 30, 'psi': 1e-5, 'c': 0.00, 'tensile_cutoff': 10000}
d_stran = -0.0001

strain = np.zeros(6)
d_strain = np.zeros(6)
#
# delta_strain = np.array([0,0,-0.0001,0,0,0])
incr_stress_input = np.array([0,0,0,0,0,0])
control_type =[1,1,0,1,1,1]

# #
stress_vector = np.array([-1,-1,-1,0,0,0])
orig_stress_vector = np.copy(stress_vector)
all_stress_vectors = [stress_vector]
all_strain_vectors = [np.zeros(6)]


delta_stress = np.zeros(6)
strain_vector = np.zeros(6)


for t in range(10):

    delta_strain = np.array([0, 0, -0.0001, 0, 0, 0])
    correction_delta_strain = np.zeros(6)

    approx_delta_stress = np.zeros(6)
    old_d_stress_vector = np.copy(stress_vector)

    max_iter = 100
    for i in range(max_iter+1):

        u_d_stress =  np.zeros(6) - approx_delta_stress
        incr_driver = IncrDriver(correction_delta_strain, u_d_stress, control_type, params)
        delta_strain = delta_strain + correction_delta_strain

        stress_tensor = StressUtils.calculate_stress_tensor_static(stress_vector)

        model = MatsuokaNakai(stress_tensor, delta_strain, params, t)
        model.main()

        if i <max_iter:
            approx_delta_stress = model.stress_utils.stress_vector - old_d_stress_vector
            stress_vector = np.copy(old_d_stress_vector)
        else:
            strain_vector = strain_vector + delta_strain
            stress_vector = model.stress_utils.stress_vector


    all_stress_vectors.append(stress_vector)
    all_strain_vectors.append(strain_vector)


a=1+1

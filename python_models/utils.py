import numpy as np

class StressUtils:
    def __init__(self):

        self.stress_vector = None
        self.stress_tensor = None
        self.deviatoric_stress_tensor = None
        self.principle_stresses = None

        # self.stress_tensor = stress_tensor
        # self.stress_vector = np.array([stress_tensor[0,0], stress_tensor[1,1], stress_tensor[2,2], stress_tensor[0,1], stress_tensor[1,2], stress_tensor[2,0]])


        self.first_invariant = None
        self.second_invariant = None
        self.third_invariant = None

        self.first_deviatoric_invariant = None
        self.second_deviatoric_invariant = None
        self.third_deviatoric_invariant = None

        self.mean_stress = None
        self.zeta_invariant = None
        self.rho_invariant = None

        self.lode_angle = None

        self.derivative_second_deviatoric_invariant = None
        self.derivative_third_deviatoric_invariant = None
        self.derivative_mean_stress = None
        self.derivative_lode_angle = None

    def update_stress(self, stress_vector):
        self.stress_vector = stress_vector
        self.stress_tensor = StressUtils.calculate_stress_tensor_static(stress_vector)


    def update_invariants(self):
        self.first_invariant = self.get_first_invariant()
        self.second_invariant = self.get_second_invariant()
        self.third_invariant = self.get_third_invariant()

        self.first_deviatoric_invariant = self.get_first_deviatoric_invariant()
        self.second_deviatoric_invariant = self.get_second_deviatoric_invariant()
        self.third_deviatoric_invariant = self.get_third_deviatoric_invariant()

        self.mean_stress = self.get_mean_stress()
        self.zeta_invariant = self.get_zeta_invariant()
        self.rho_invariant = self.get_rho_invariant()

        self.lode_angle = self.get_lode_angle()

    def update_derivatives(self):
        self.get_derivative_second_deviatoric_invariant()
        self.get_derivative_third_deviatoric_invariant()
        self.get_derivative_mean_stress()
        self.get_derivative_lode_angle()

    # def calculate_stress_tensor(self):
    #     self.stress_tensor = np.array([[self.stress_vector[0], self.stress_vector[3], self.stress_vector[5]],
    #                                     [self.stress_vector[3], self.stress_vector[1], self.stress_vector[4]],
    #                                     [self.stress_vector[5], self.stress_vector[4], self.stress_vector[2]]])
    @staticmethod
    def calculate_stress_tensor_static(stress_vector):
        return np.array([[stress_vector[0], stress_vector[3], stress_vector[5]],
                         [stress_vector[3], stress_vector[1], stress_vector[4]],
                         [stress_vector[5], stress_vector[4], stress_vector[2]]])

    @staticmethod
    def calculate_stress_vector_static(stress_tensor):
        return np.array([stress_tensor[0,0], stress_tensor[1,1], stress_tensor[2,2], stress_tensor[0,1], stress_tensor[1,2], stress_tensor[2,0]])


    def calculate_principle_stresses(self):

        # calculate principle stresses
        self.principle_stresses = np.linalg.eigvals(self.stress_tensor)

    def calculate_principle_stresses_invariants(self):
        sigma1 = self.mean_stress + np.sqrt(4*self.second_deviatoric_invariant/3) * np.cos(self.lode_angle)
        sigma2 = self.mean_stress + np.sqrt(4*self.second_deviatoric_invariant/3) * np.cos(self.lode_angle - 2*np.pi/3)
        sigma3 = self.mean_stress + np.sqrt(4*self.second_deviatoric_invariant/3) * np.cos(self.lode_angle + 2*np.pi/3)


        self.principle_stresses = np.array([sigma1, sigma2, sigma3])

    # fundamental invariants
    def get_first_invariant(self):
        return np.trace(self.stress_tensor)

    def get_first_invariant_principal_stresses(self):
        return np.sum(self.principle_stresses)

    def get_second_invariant(self):
        return 0.5 * (self.first_invariant**2 - np.trace(self.stress_tensor @ self.stress_tensor))

    def get_second_invariant_principal_stresses(self):
        return self.principle_stresses[0] * self.principle_stresses[1] + self.principle_stresses[1] * self.principle_stresses[2] + self.principle_stresses[2] * self.principle_stresses[0]

    def get_third_invariant(self):
        return np.linalg.det(self.stress_tensor)

    def get_third_invariant_principal_stresses(self):
        return self.principle_stresses[0] * self.principle_stresses[1] * self.principle_stresses[2]

    # deviatoric invariants
    def get_deviatoric_tensor(self):
        return self.stress_tensor - np.eye(3) * self.first_invariant / 3

    def get_first_deviatoric_invariant(self):
        return 0

    def get_second_deviatoric_invariant(self):
        self.deviatoric_stress_tensor =  self.get_deviatoric_tensor()
        # return 0.5 * (np.trace(self.get_deviatoric_tensor())**2 - np.trace(self.get_deviatoric_tensor() @ self.get_deviatoric_tensor()))
        return 0.5 * np.sum(np.square(self.deviatoric_stress_tensor))

    def get_derivative_second_deviatoric_invariant(self):
        """
        Equal to deviatoric stress vector
        :return:
        """
        dj2dS =np.zeros(6)

        dj2dS[0] = self.stress_vector[0] - self.mean_stress
        dj2dS[1] = self.stress_vector[1] - self.mean_stress
        dj2dS[2] = self.stress_vector[2] - self.mean_stress
        dj2dS[3] = self.stress_vector[3] * 2
        dj2dS[4] = self.stress_vector[4] * 2
        dj2dS[5] = self.stress_vector[5] * 2

        self.derivative_second_deviatoric_invariant = dj2dS


    def get_third_deviatoric_invariant(self):
        return np.linalg.det(self.get_deviatoric_tensor())

    def get_derivative_third_deviatoric_invariant(self):
        j2 = self.second_deviatoric_invariant
        dj3dS = np.zeros(6)
        dj3dS[0] = (self.stress_vector[0] - self.mean_stress) ** 2 + self.stress_vector[5]**2 + self.stress_vector[3]**2 - 2* j2/3
        dj3dS[1] = (self.stress_vector[1] - self.mean_stress) ** 2 + self.stress_vector[3]**2 + self.stress_vector[4]**2 - 2* j2/3
        dj3dS[2] = (self.stress_vector[2] - self.mean_stress) ** 2 + self.stress_vector[4]**2 + self.stress_vector[5]**2 - 2* j2/3
        dj3dS[3] = 2 * ((self.stress_vector[0] + self.stress_vector[1] - 2*self.mean_stress) * self.stress_vector[3] + self.stress_vector[5] * self.stress_vector[4])
        dj3dS[4] = 2 * ((self.stress_vector[1] + self.stress_vector[2] - 2*self.mean_stress) * self.stress_vector[4] + self.stress_vector[3] * self.stress_vector[5])
        dj3dS[5] = 2 * ((self.stress_vector[2] + self.stress_vector[0] - 2*self.mean_stress) * self.stress_vector[5] + self.stress_vector[3] * self.stress_vector[4])

        self.derivative_third_deviatoric_invariant = dj3dS

    def get_mean_stress(self):
        return self.first_invariant / 3

    def get_derivative_mean_stress(self):

        self.derivative_mean_stress = np.array([1/3, 1/3, 1/3, 0,0,0])

    def get_zeta_invariant(self):
        return np.sqrt(3) *self.mean_stress

    def get_rho_invariant(self):
        return np.sqrt(2) * np.sqrt(self.second_deviatoric_invariant)

    def get_lode_angle(self):

        J2 = self.second_deviatoric_invariant
        J3 = self.third_deviatoric_invariant

        invariant = np.sqrt(27) * J3 / (2 * (J2**(3/2)))

        # check if invariant is close to 1
        if np.isclose(invariant,1.0):
            invariant = 1
        elif np.isclose(invariant,-1.0):
            invariant = -1

        return 1/3 * np.arcsin(invariant)

    def __get_derivative_lode_angle_J2(self):
        J2 = self.second_deviatoric_invariant
        J3 = self.third_deviatoric_invariant


        dtheta_dJ2 = -3 ** (3 / 2) * J3 / (4 * J2**(5/2) * np.sqrt(1 - (27 * J3 ** 2/(4*J2**3))))


        return dtheta_dJ2

    def __get_derivative_lode_angle_J3(self):
        J2 = self.second_deviatoric_invariant
        J3 = self.third_deviatoric_invariant
        dtheta_dJ3 = np.sqrt(3) / (np.sqrt(4*J2**3 -27*J3**2))

        return dtheta_dJ3

    def get_derivative_lode_angle(self):

        if np.isclose(self.lode_angle,np.pi/6) or np.isclose(self.lode_angle,-np.pi/6):
            self.derivative_lode_angle = np.zeros(6)
            return self.derivative_lode_angle

        dtheta_dsigma = (self.__get_derivative_lode_angle_J2() * self.get_derivative_second_deviatoric_invariant()
                         + self.__get_derivative_lode_angle_J3() * self.get_derivative_third_deviatoric_invariant())

        if np.isnan(dtheta_dsigma).any():
            raise ValueError("Nan values in derivative of lode angle")


        self.derivative_lode_angle = dtheta_dsigma


    # octahedral plane / pi-plane
    def get_octahedral_stress_normal(self):
        return self.first_invariant / 3

    def get_octahedral_stress_shear(self):
        return np.sqrt(2/3 * self.second_deviatoric_invariant)



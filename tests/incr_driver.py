import numpy as np

from tests.utils import Utils


VERTICAL_AXIS_INDEX = 1 # [0, 1, 2] = [x, y, z]
NDIM = 3

class IncrDriver:
    def __init__(self,initial_stress, strain_increment, stress_increment, constitutive_model_info, n_time_steps,
                 max_iterations, voigt_size=6):
        """
        Initialize the IncrDriver class.

        Args:
            initial_stress (np.ndarray): The initial stress vector.
            strain_increment (np.ndarray): The strain increment vector.
            stress_increment (np.ndarray): The stress increment vector.
            constitutive_model_info (dict): Information about the constitutive model, including language and file name.
            n_time_steps (int): The number of time steps to solve.
            max_iterations (int): The maximum number of iterations per time step for the solver.
            voigt_size (int): The size of the Voigt notation for stress/strain vectors. Default is 6 (3D).

        """

        self.initial_stress = initial_stress
        self.strain_increment = strain_increment
        self.stress_increment = stress_increment
        self.constitutive_model_info = constitutive_model_info
        self.n_time_steps = n_time_steps
        self.max_iterations = max_iterations
        self.voigt_size = voigt_size

        self.stresses = []
        self.strains = []

    def solve_oedometer_strain_controlled(self):
        """
        Solve the oedometer problem with strain control.

        """
        control_type = np.zeros(self.voigt_size)

        self.solve(control_type)

    def solve_oedometer_stress_controlled(self):
        """
        Solve the oedometer problem with stress control.

        vertical axis is controlled

        """

        control_type = np.zeros(self.voigt_size)
        control_type[VERTICAL_AXIS_INDEX] = 1

        self.solve(control_type)

    def solve_triaxial_strain_controlled(self):
        """
        Solve the triaxial problem with strain control.

        vertical axis is controlled

        """

        control_type = np.zeros(self.voigt_size)

        control_type[:NDIM] = 1
        control_type[VERTICAL_AXIS_INDEX] = 0

        self.solve(control_type)

    def solve_triaxial_stress_controlled(self):
        """
        Solve the triaxial problem with stress control.

        """

        control_type = np.zeros(self.voigt_size)
        control_type[:NDIM] = 1

        self.solve(control_type)


    def solve(self, control_type):
        """
        Solve the incremental problem using the specified constitutive model.
        The method iteratively updates the stress and strain vectors based on the provided strain increment and
        control type.

        control_type (list): A list indicating the control type for each component (1 for stress, 0 for strain).
        """



        language = self.constitutive_model_info['language']
        stress_updated, ddsdde, statev_updated = None, None, None

        state_variables = self.constitutive_model_info['state_vars']

        # run umat in order to retrieve the elastic matrix
        if language == "c":
            _, ddsdde, _ = Utils.run_c_umat(self.constitutive_model_info['file_name'],
                                            self.initial_stress,
                                            np.copy(state_variables),
                                            np.zeros(self.voigt_size),
                                            np.zeros(self.voigt_size),
                                            self.constitutive_model_info["properties"], 0)


        elif language == "fortran":
            _, ddsdde, _ = Utils.run_fortran_umat(self.constitutive_model_info['file_name'],
                                                  self.initial_stress,
                                                  np.copy(state_variables),
                                                  np.zeros(self.voigt_size),
                                                  np.zeros(self.voigt_size),
                                                  self.constitutive_model_info["properties"], 0)

        else:
            ValueError(f"Language {language} not supported. Only 'c' and 'fortran' are supported.")

        # initialize stress and strain vectors
        strain_vector = np.zeros(self.voigt_size)
        stress_vector = np.copy(self.initial_stress)

        stresses = []
        strains = []

        # loop over time steps
        for t in range(self.n_time_steps):

            delta_strain = np.copy(self.strain_increment)
            correction_delta_strain = np.zeros(self.voigt_size)
            approx_delta_stress = np.zeros(self.voigt_size)

            # save the stress vector of the previous time step
            old_d_stress_vector = np.copy(stress_vector)

            # save the state variable vector of the previous time step
            prev_state_variables = np.copy(state_variables)

            # loop over the maximum number of non-linear iterations
            for i in range(self.max_iterations + 1):

                # get stress increment at stress controlled components, else zero
                stress_increment = np.where(control_type, self.stress_increment, 0.0)

                # calculate the undesired stress
                u_d_stress = stress_increment + approx_delta_stress

                # update the delta strain using the elastic matrix and the undesired stress
                self.calculate_iteration(correction_delta_strain, u_d_stress, control_type, ddsdde)

                # correct delta strain
                delta_strain = delta_strain - correction_delta_strain

                # run constitutive model
                if language == "c":
                    stress_updated, ddsdde, state_variables = Utils.run_c_umat(self.constitutive_model_info['file_name'],
                                                                              stress_vector,
                                                                              state_variables,
                                                                              strain_vector,
                                                                              delta_strain,
                                                                              self.constitutive_model_info[
                                                                                  "properties"], t)
                elif language == "fortran":
                    stress_updated, ddsdde, state_variables = Utils.run_fortran_umat(self.constitutive_model_info['file_name'],
                                                                              stress_vector,
                                                                              state_variables,
                                                                              strain_vector,
                                                                              delta_strain,
                                                                              self.constitutive_model_info[
                                                                                  "properties"], t)
                else:
                    ValueError(f"Language {language} not supported. Only 'c' and 'fortran' are supported.")


                # if not at max iterations, update the approximation for delta stress and reset stress vector and state
                # variables
                # else update strain and stress and go to the next time step
                if i < self.max_iterations:

                    approx_delta_stress = stress_updated - old_d_stress_vector

                    # reset stress vector and state variables
                    stress_vector = np.copy(old_d_stress_vector)
                    state_variables = np.copy(prev_state_variables)
                else:
                    strain_vector = strain_vector + delta_strain
                    stress_vector = stress_updated

            stresses.append(stress_vector)
            self.stresses = np.array(stresses)

            strains.append(strain_vector)
            self.strains = np.array(strains)

    @staticmethod
    def calculate_iteration(strains, stresses, control_type, elastic_matrix):
        """
        Solves the system of equations to update the strain and stress vectors based on the control type and elastic
        matrix.

        Args:
            strains (np.ndarray): The strain vector to be updated.
            stresses (np.ndarray): The stress vector to be updated.
            control_type (list): A list indicating the control type for each component (1 for stress, 0 for strain).
            elastic_matrix (np.ndarray): The elastic matrix used for calculations.
        """
        D = elastic_matrix

        # prevent NANs, if the matrix is singular, fill the diagonal with a small number
        diag_indices = np.diag_indices_from(D)
        small_diagonal = np.abs(D[diag_indices]) < 1e-12
        D[diag_indices][small_diagonal] = 1e-12  # or use np.finfo(float).eps for machine epsilon

        control_type_array = np.array(control_type)
        strain_controlled = control_type_array == 0
        stress_controlled = control_type_array == 1

        # Apply strain-controlled updates to stress
        if np.any(strain_controlled):
            delta_stress = D[:, strain_controlled] @ strains[strain_controlled]
            stresses_copy = stresses - delta_stress
        else:
            stresses_copy = stresses.copy()

        # if stresses are controlled, calculate strains at stress controlled indices,
        # i.e. solve eps = D^-1 * sigma
        if np.any(stress_controlled):
            reduced_D = D[np.ix_(stress_controlled, stress_controlled)]
            reduced_stresses = stresses_copy[stress_controlled]

            strains[stress_controlled] = np.linalg.solve(reduced_D, reduced_stresses)

        # Update full stresses with updated strains
        if np.any(strain_controlled):
            stresses[strain_controlled] = D[strain_controlled, :] @ strains

import numpy as np

class IncrDriver():
    def __init__(self,initial_stress, strain_increment, control_type, constitutive_model, n_time_steps,max_iterations, voight_size=6):
        self.initial_stress = initial_stress
        self.strain_increment = strain_increment
        self.control_type = control_type
        self.constitutive_model = constitutive_model
        self.n_time_steps = n_time_steps
        self.max_iterations = max_iterations

        self.voight_size = voight_size

        self.stresses = []
        self.strains = []


    def solve(self):

        strain_vector = np.zeros(self.voight_size)
        stress_vector = np.copy(self.initial_stress)

        self.constitutive_model.elasticity_model.calculate_elastic_matrix()


        for t in range(self.n_time_steps+1):
            delta_strain = np.copy(self.strain_increment)
            correction_delta_strain = np.zeros(self.voight_size)
            approx_delta_stress = np.zeros(self.voight_size)
            old_d_stress_vector = np.copy(stress_vector)

            for i in range(self.max_iterations + 1):
                #undesired stress
                u_d_stress = np.zeros(self.voight_size) - approx_delta_stress
                self.calculate_iteration(correction_delta_strain, u_d_stress, self.control_type, self.constitutive_model.elasticity_model.elastic_matrix)

                # correct delta strain
                delta_strain = delta_strain + correction_delta_strain

                # run constitutive model
                self.constitutive_model.main(stress_vector, delta_strain)

                # if not at max iterations, update the approximation for delta stress
                # else update strain and stress and go to the next time step
                if i < self.max_iterations:
                    approx_delta_stress = self.constitutive_model.stress_utils.stress_vector - old_d_stress_vector
                    stress_vector = np.copy(old_d_stress_vector)
                else:
                    strain_vector = strain_vector + delta_strain
                    stress_vector = self.constitutive_model.stress_utils.stress_vector
            self.stresses.append(stress_vector)
            self.strains.append(strain_vector)


    def calculate_iteration(self, strains, stresses, control_type, elastic_matrix):


        D = elastic_matrix

        stresses_copy = stresses.copy()
        # strain control_type
        for i in range(len(control_type)):
            if control_type[i] ==0:
                stresses_copy = stresses_copy - strains[i] * D[:,i]

        indices = np.nonzero(np.array(control_type)==1)[0]
        reduced_D = D[np.ix_(indices,indices)]
        reduced_stresses = stresses_copy[indices]

        reduced_strains = np.linalg.solve(reduced_D,reduced_stresses)

        strains[indices] = reduced_strains

        for i in range(len(control_type)):
            if control_type[i] ==1:
                stresses[i] = D[i,:].dot(strains)








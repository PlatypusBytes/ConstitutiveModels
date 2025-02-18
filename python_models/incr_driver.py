import numpy as np

class IncrDriver():
    def __init__(self, ):



        pass


    def calculate_time_step(self, strains, stresses, control_type, elastic_matrix):


        D = elastic_matrix
        # self.D_inv = np.linalg.inv(self.D)

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








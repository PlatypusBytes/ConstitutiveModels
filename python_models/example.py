import numpy as np

from python_models.incr_driver import IncrDriver
from python_models.yield_surfaces.matsuoka_nakai import MatsuokaNakai
from python_models.elasticity_models.hookes_law import HooksLaw
from python_models.constitutive_model import BaseConstitutiveModel
from python_models.incr_driver import IncrDriver

yield_function = MatsuokaNakai({'angle': 30, 'cohesion': 0.00})
flow_rule = MatsuokaNakai({'angle': 0, 'cohesion': 0.00})

G = 2000
nu = 0.33
E = 2 * G * (1 + nu)
elasticity_model = HooksLaw({'young_modulus': E, 'poison_ratio': nu})

stress_vector = np.array([-1, -1, -1, 0, 0, 0])
strain_increment =[0,0,-0.0001,0,0,0]

incr_stress_input = np.array([0, 0, 0, 0, 0, 0])
control_type = [1, 1, 0, 1, 1, 1]

constitutive_model = BaseConstitutiveModel(yield_function, flow_rule, elasticity_model)

incremental_driver = IncrDriver(stress_vector, strain_increment, control_type, constitutive_model, 10, 100)
incremental_driver.solve()

a=1+1

from abc import abstractmethod

class BaseLaw():

    def __init__(self):
        pass

    @abstractmethod
    def calculate_elastic_matrix(self):
        raise NotImplementedError("Subclass must implement abstract method 'calculate_elastic_matrix'")


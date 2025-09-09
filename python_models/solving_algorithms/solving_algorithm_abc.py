
from abc import ABC, abstractmethod

class SolvingAlgorithmABC():

    def __init__(self):
        pass

    @abstractmethod
    def solve(self) :
        raise NotImplementedError("Subclass must implement abstract method 'solve'")

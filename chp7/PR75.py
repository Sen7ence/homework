import numpy as np


class PR75:
    def __init__(self, Tc1, pc1, omega1, Tc2, pc2, omega2, kij):
        self.Tc1 = Tc1
        self.pc1 = pc1
        self.omega1 = omega1
        self.Tc2 = Tc2
        self.pc2 = pc2
        self.omega2 = omega2
        self.kij = kij

    R = 8.314462618  # J/(mol·K)

import numpy as np


class PRHS2:
    def __init__(self, Tc1, pc1, omega1, M1, x1, Tc2, pc2, omega2, M2, x2, kij, ps0):
        self.Tc1 = Tc1  # K
        self.pc1 = pc1 * 1e6  # Pa，输入MPa
        self.omega1 = omega1  # 无量纲
        self.M1 = M1 / 1e3  # kg/mol，输入g/mol
        self.x1 = x1  # 组分1的摩尔分数

        self.Tc2 = Tc2  # K
        self.pc2 = pc2 * 1e6  # Pa，输入MPa
        self.omega2 = omega2  # 无量纲
        self.M2 = M2 / 1e3  # kg/mol，输入g/mol
        self.x2 = x2  # 组分2的摩尔分数

        self.ps0 = ps0  # MPa

        self.kij = kij  # 无量纲

    R = 8.314462618  # J/(mol*K)

    # 计算a和b

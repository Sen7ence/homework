import numpy as np


class PR415:
    def __init__(self, Tc1, pc1, omega1, M1, x1, Tc2, pc2, omega2, M2, kij, ps0):
        self.Tc1 = Tc1  # K
        self.pc1 = pc1 * 1e6  # Pa，输入MPa
        self.omega1 = omega1  # 无量纲
        self.M1 = M1 / 1e3  # kg/mol，输入g/mol
        self.x1 = x1  # 组分1的摩尔分数

        self.Tc2 = Tc2  # K
        self.pc2 = pc2 * 1e6  # Pa，输入MPa
        self.omega2 = omega2  # 无量纲
        self.M2 = M2 / 1e3  # kg/mol，输入g/mol
        self.x2 = 1 - x1  # 组分2的摩尔分数

        self.ps0 = ps0  # MPa

        self.kij = kij  # 无量纲

    R = 8.314462618  # J/(mol*K)

    # 计算a和b
    def params(self, T):
        kappa1 = 0.37464 + 1.54226 * self.omega1 - 0.26992 * self.omega1**2
        kappa2 = 0.37464 + 1.54226 * self.omega2 - 0.26992 * self.omega2**2
        Tr1 = T / self.Tc1
        Tr2 = T / self.Tc2
        alpha1 = (1 + kappa1 * (1 - Tr1**0.5)) ** 2
        alpha2 = (1 + kappa2 * (1 - Tr2**0.5)) ** 2
        a1 = 0.45724 * self.R**2 * self.Tc1**2 / self.pc1 * alpha1
        a2 = 0.45724 * self.R**2 * self.Tc2**2 / self.pc2 * alpha2
        da1 = (
            -0.45724
            * self.R**2
            * self.Tc1**2
            / self.pc1
            * kappa1
            * (1 + kappa1 * (1 - Tr1**0.5))
            * (Tr1**-0.5)
            / self.Tc1
        )
        da2 = (
            -0.45724
            * self.R**2
            * self.Tc2**2
            / self.pc2
            * kappa2
            * (1 + kappa2 * (1 - Tr2**0.5))
            * (Tr2**-0.5)
            / self.Tc2
        )
        b1 = 0.07780 * self.R * self.Tc1 / self.pc1
        b2 = 0.07780 * self.R * self.Tc2 / self.pc2

        a = (
            self.x1**2 * a1
            + self.x2**2 * a2
            + 2 * self.x1 * self.x2 * (a1 * a2) ** 0.5 * (1 - self.kij)
        )
        b = self.x1 * b1 + self.x2 * b2
        da = (
            self.x1**2 * da1
            + self.x2**2 * da2
            + self.x1
            * self.x2
            * (1 - self.kij)
            * ((a2 / a1) ** 0.5 * da1 + (a1 / a2) ** 0.5 * da2)
        )
        return a, b, da

    # 计算A和B
    def AB(self, T, p):
        a, b, da = self.params(T)
        A = a * p * 1e6 / (self.R * T) ** 2  # p从MPa转换为Pa
        B = b * p * 1e6 / (self.R * T)  # p从MPa转换为Pa
        return A, B

    # 计算C2, C1, C0
    def C(self, T, p):
        A, B = self.AB(T, p)
        C2 = -(1 - B)
        C1 = A - 3 * B**2 - 2 * B
        C0 = -(A * B - B**2 - B**3)
        return C2, C1, C0

    # 计算压缩因子Z
    # 液相
    def Zl(self, T, p):
        C2, C1, C0 = self.C(T, p)
        # 牛顿法求解Z
        Zl = 0.001  # 初始猜测值
        for _ in range(100):
            f = Zl**3 + C2 * Zl**2 + C1 * Zl + C0
            df = 3 * Zl**2 + 2 * C2 * Zl + C1
            Zl_new = Zl - f / df
            if abs(Zl_new - Zl) < 1e-6:
                break
            Zl = Zl_new
        return Zl

    # 气相
    def Zg(self, T, p):
        C2, C1, C0 = self.C(T, p)
        # 牛顿法求解Z
        Zg = 1.1  # 初始猜测值
        for _ in range(100):
            f = Zg**3 + C2 * Zg**2 + C1 * Zg + C0
            df = 3 * Zg**2 + 2 * C2 * Zg + C1
            Zg_new = Zg - f / df
            if abs(Zg_new - Zg) < 1e-6:
                break
            Zg = Zg_new
        return Zg

    # 计算比体积v
    # 液相
    def vl(self, T, p):
        Zl = self.Zl(T, p)
        vl = Zl * self.R * T / (p * 1e6)  # p从MPa转换为Pa
        return vl

    # 气相
    def vg(self, T, p):
        Zg = self.Zg(T, p)
        vg = Zg * self.R * T / (p * 1e6)  # p从MPa转换为Pa
        return vg

    # 计算焓的余函数
    # 液相
    def h_res_l(self, T, p):
        a, b, da = self.params(T)
        Zl = self.Zl(T, p)
        vl = self.vl(T, p)
        hr_l = (T * da - a) / (b * np.sqrt(8)) * np.log(
            (vl - 0.414 * b) / (vl + 2.414 * b)
        ) + self.R * T * (1 - Zl)
        return hr_l

    # 气相
    def h_res_g(self, T, p):
        a, b, da = self.params(T)
        Zg = self.Zg(T, p)
        vg = self.vg(T, p)
        hr_g = (T * da - a) / (b * np.sqrt(8)) * np.log(
            (vg - 0.414 * b) / (vg + 2.414 * b)
        ) + self.R * T * (1 - Zg)
        return hr_g

    # 计算熵的余函数
    # 液相
    def s_res_l(self, T, p):
        a, b, da = self.params(T)
        vl = self.vl(T, p)
        sr_l = (
            -self.R * np.log((vl - b) / vl)
            - self.R * np.log(vl / (self.R * T / (p * 1e6)))
            + da / (b * np.sqrt(8)) * np.log((vl - 0.414 * b) / (vl + 2.414 * b))
        )
        return sr_l

    # 气相
    def s_res_g(self, T, p):
        a, b, da = self.params(T)
        vg = self.vg(T, p)
        sr_g = (
            -self.R * np.log((vg - b) / vg)
            - self.R * np.log(vg / (self.R * T / (p * 1e6)))
            + da / (b * np.sqrt(8)) * np.log((vg - 0.414 * b) / (vg + 2.414 * b))
        )
        return sr_g

    # 计算c_p积分
    def cp(self, T, A1, A2, B1, B2, C1, C2, D1, D2):
        cp = (
            (A1 + A2) * 0.5 * (T - 273.15)
            + (B1 + B2) * 0.5 / 2 * (T**2 - 273.15**2)
            + (C1 + C2) * 0.5 / 3 * (T**3 - 273.15**3)
            + (D1 + D2) * 0.5 / 4 * (T**4 - 273.15**4)
        )
        return cp

    # 计算c_p/T积分
    def cpT(self, T, A1, A2, B1, B2, C1, C2, D1, D2):
        cp = (
            (A1 + A2) * 0.5 * np.log(T / 273.15)
            + (B1 + B2) * 0.5 * (T - 273.15)
            + (C1 + C2) * 0.5 / 2 * (T**2 - 273.15**2)
            + (D1 + D2) * 0.5 / 3 * (T**3 - 273.15**3)
        )
        return cp

    # 计算焓和熵
    # 液相
    def h_l(self, T, A1, A2, B1, B2, C1, C2, D1, D2, p):
        h_r_ps_0 = self.h_res_l(273.15, self.ps0)
        cp0 = self.cp(T, A1, A2, B1, B2, C1, C2, D1, D2)
        h_res_l = self.h_res_l(T, p)
        hl = (
            200 * 1e3
            + cp0
            + (h_r_ps_0 - h_res_l) / (self.x1 * self.M1 + self.x2 * self.M2)
        )  # J/kg
        return hl

    def s_l(self, T, A1, A2, B1, B2, C1, C2, D1, D2, p):
        s_r_ps_0 = self.s_res_l(273.15, self.ps0)
        cpT = self.cpT(T, A1, A2, B1, B2, C1, C2, D1, D2)
        sr_l = self.s_res_l(T, p)
        sl = (
            1e3
            + cpT
            + (s_r_ps_0 - self.R * np.log(p / self.ps0) - sr_l)
            / (self.x1 * self.M1 + self.x2 * self.M2)
        )  # J/(kg*K)
        return sl

    # 气相
    def h_g(self, T, A1, A2, B1, B2, C1, C2, D1, D2, p):
        h_r_ps_0 = self.h_res_l(273.15, self.ps0)  # 使用液相作为基准
        cp0 = self.cp(T, A1, A2, B1, B2, C1, C2, D1, D2)
        h_res_g = self.h_res_g(T, p)
        hg = (
            200 * 1e3
            + cp0
            + (h_r_ps_0 - h_res_g) / (self.x1 * self.M1 + self.x2 * self.M2)
        )  # J/kg
        return hg

    def s_g(self, T, A1, A2, B1, B2, C1, C2, D1, D2, p):
        s_r_ps_0 = self.s_res_l(273.15, self.ps0)  # 使用液相作为基准
        cpT = self.cpT(T, A1, A2, B1, B2, C1, C2, D1, D2)
        sr_g = self.s_res_g(T, p)
        sg = (
            1e3
            + cpT
            + (s_r_ps_0 - self.R * np.log(p / self.ps0) - sr_g)
            / (self.x1 * self.M1 + self.x2 * self.M2)
        )  # J/(kg*K)
        return sg


R290R600a = PR415(
    Tc1=369.89,
    pc1=4.2512,
    omega1=0.1521,
    M1=44.096,  # R290
    x1=0.5,
    Tc2=407.81,
    pc2=3.629,
    omega2=0.184,
    M2=58.122,  # R600a
    kij=0.064,
    ps0=0.32979,
)

# 300K下计算比焓和比熵
print(
    R290R600a.h_l(
        300,
        -95.80,
        -23.91,
        6.945,
        6.605,
        -3.597 * 1e-3,
        -3.176 * 1e-3,
        7.290 * 1e-7,
        4.981 * 1e-7,
        1.0,
    )
)
print(
    R290R600a.s_l(
        300,
        -95.80,
        -23.91,
        6.945,
        6.605,
        -3.597 * 1e-3,
        -3.176 * 1e-3,
        7.290 * 1e-7,
        4.981 * 1e-7,
        1.0,
    )
)

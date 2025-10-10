class FRFluid2:
    def __init__(self, Tc1, Tc2, pc1, pc2, omega1, omega2, M1, M2, x1, kij):
        self.Tc1 = Tc1  # K
        self.Tc2 = Tc2  # K
        self.pc1 = pc1 * 1e6  # Pa，输入MPa，改为乘法
        self.pc2 = pc2 * 1e6  # Pa，输入MPa，改为乘法
        self.omega1 = omega1  # 无量纲
        self.omega2 = omega2  # 无量纲
        self.M1 = M1 / 1e3  # kg/mol，输入g/mol
        self.M2 = M2 / 1e3  # kg/mol，输入g/mol
        self.x1 = x1  # 组分1的摩尔分率
        self.x2 = 1 - x1  # 组分2的摩尔分率
        self.kij = kij  # 组分间的二元交互作用参数

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
        b1 = 0.07780 * self.R * self.Tc1 / self.pc1
        b2 = 0.07780 * self.R * self.Tc2 / self.pc2
        a = (
            self.x1**2 * a1
            + self.x2**2 * a2
            + 2 * self.x1 * self.x2 * (a1 * a2) ** 0.5 * (1 - self.kij)
        )
        b = self.x1 * b1 + self.x2 * b2
        return a, b

    # 计算A和B
    def AB(self, T, p):
        a, b = self.params(T)
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
        Zg = 1.0  # 初始猜测值
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
        vl = (
            Zl * self.R * T / (p * 1e6 * (self.x1 * self.M1 + self.x2 * self.M2))
        )  # p从MPa转换为Pa
        return vl

    # 气相
    def vg(self, T, p):
        Zg = self.Zg(T, p)
        vg = (
            Zg * self.R * T / (p * 1e6 * (self.x1 * self.M1 + self.x2 * self.M2))
        )  # p从MPa转换为Pa
        return vg


R290_R600a_1 = FRFluid2(
    369.89, 407.81, 4.2512, 3.629, 0.1521, 0.184, 44.096, 58.122, 0.5, 0.064
)
R290_R600a_2 = FRFluid2(
    369.89, 407.81, 4.2512, 3.629, 0.1521, 0.184, 44.096, 58.122, 0.5, 0.1
)
R290_R600a_3 = FRFluid2(
    369.89, 407.81, 4.2512, 3.629, 0.1521, 0.184, 44.096, 58.122, 0.5, 0
)
R290_R600a_4 = FRFluid2(
    369.89, 407.81, 4.2512, 3.629, 0.1521, 0.184, 44.096, 58.122, 0.5, -0.1
)
print(R290_R600a_1.vg(300, 0.1))
print(R290_R600a_2.vg(300, 0.1))
print(R290_R600a_3.vg(300, 0.1))
print(R290_R600a_4.vg(300, 0.1))
print(R290_R600a_1.vg(300, 0.2))
print(R290_R600a_2.vg(300, 0.2))
print(R290_R600a_3.vg(300, 0.2))
print(R290_R600a_4.vg(300, 0.2))
print(R290_R600a_1.vg(300, 0.3))
print(R290_R600a_2.vg(300, 0.3))
print(R290_R600a_3.vg(300, 0.3))
print(R290_R600a_4.vg(300, 0.3))

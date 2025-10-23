import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os

# 使用 Times New Roman 作为 matplotlib 全局字体
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14  # 增大全局字体
plt.rcParams["axes.linewidth"] = 1.5  # 增粗坐标轴


class PR74:
    def __init__(self, Tc, pc, omega):
        self.Tc = Tc
        self.pc = pc * 1e6
        self.omega = omega

    R = 8.314462618  # J/(mol*K)

    def params(self, T):
        kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
        Tr = T / self.Tc
        alpha = (1 + kappa * (1 - Tr**0.5)) ** 2
        a = 0.45724 * self.R**2 * self.Tc**2 / self.pc * alpha
        b = 0.07780 * self.R * self.Tc / self.pc
        da = (
            -0.45724
            * self.R**2
            * self.Tc**2
            / self.pc
            * kappa
            * (1 + kappa * (1 - Tr**0.5))
            * (Tr**-0.5)
            / self.Tc
        )
        return a, b, da

    def AB(self, T, p):
        a, b, da = self.params(T)
        A = a * p * 1e6 / (self.R * T) ** 2
        B = b * p * 1e6 / (self.R * T)
        return A, B

    def C(self, T, p):
        A, B = self.AB(T, p)
        C2 = -(1 - B)
        C1 = A - 3 * B**2 - 2 * B
        C0 = -(A * B - B**2 - B**3)
        return C2, C1, C0

    def Zl(self, T, p):
        C2, C1, C0 = self.C(T, p)
        Zl = 0.001
        for _ in range(100):
            f = Zl**3 + C2 * Zl**2 + C1 * Zl + C0
            df = 3 * Zl**2 + 2 * C2 * Zl + C1
            Zl_new = Zl - f / df
            if abs(Zl_new - Zl) < 1e-6:
                break
            Zl = Zl_new
        return Zl

    def Zg(self, T, p):
        C2, C1, C0 = self.C(T, p)
        Zg = 1.1
        for _ in range(100):
            f = Zg**3 + C2 * Zg**2 + C1 * Zg + C0
            df = 3 * Zg**2 + 2 * C2 * Zg + C1
            Zg_new = Zg - f / df
            if abs(Zg_new - Zg) < 1e-6:
                break
            Zg = Zg_new
        return Zg

    def vl(self, T, p):
        Zl = self.Zl(T, p)
        vl = Zl * self.R * T / (p * 1e6)
        return vl

    def vg(self, T, p):
        Zg = self.Zg(T, p)
        vg = Zg * self.R * T / (p * 1e6)
        return vg

    def phi_l(self, T, p):
        Zl = self.Zl(T, p)
        a, b, da = self.params(T)
        A, B = self.AB(T, p)
        ln_phi_l = (
            Zl
            - 1
            - np.log(Zl - B)
            - A
            / (2 * (2**0.5) * B)
            * np.log((Zl + (1 + 2**0.5) * B) / (Zl + (1 - 2**0.5) * B))
        )
        phi_l = np.exp(ln_phi_l)
        return phi_l

    def phi_g(self, T, p):
        Zg = self.Zg(T, p)
        a, b, da = self.params(T)
        A, B = self.AB(T, p)
        ln_phi_g = (
            Zg
            - 1
            - np.log(Zg - B)
            - A
            / (2 * (2**0.5) * B)
            * np.log((Zg + (1 + 2**0.5) * B) / (Zg + (1 - 2**0.5) * B))
        )
        phi_g = np.exp(ln_phi_g)
        return phi_g

    # 计算饱和压力
    def psat(self, T):
        from scipy.optimize import fsolve

        def objective(p):
            return self.phi_l(T, p) - self.phi_g(T, p)

        p_initial = 0.0001  # 初始猜测值，单位MPa
        psat_solution = fsolve(objective, p_initial)
        return psat_solution[0]  # 返回平衡压力，单位MPa

    # 绘制p-T相图
    def plot_pT(self, T_min, T_max, savepath):
        T_grid = np.linspace(T_min, T_max, 100)
        p_grid = np.zeros_like(T_grid)
        for i, T in enumerate(T_grid):
            p_grid[i] = self.psat(T)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(T_grid, p_grid, "r-", linewidth=2)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.set_xlabel("T (K)")
        ax.set_ylabel("p (MPa)")
        ax.grid(True)

        base_dir = os.path.dirname(os.path.abspath(__file__))
        fig_dir = os.path.join(base_dir, "figs")
        os.makedirs(fig_dir, exist_ok=True)
        savepath = os.path.join(fig_dir, savepath)

        fig.savefig(savepath, dpi=600, bbox_inches="tight", transparent=False)
        plt.close(fig)


R290 = PR74(Tc=366.8, pc=4.248, omega=0.152)
R290.plot_pT(220, 325, "R290_pT.png")
R600a = PR74(Tc=407.8, pc=3.796, omega=0.227)
R600a.plot_pT(220, 325, "R600a_pT.png")
R1234yf = PR74(Tc=367.85, pc=3.3822, omega=0.276)
R1234yf.plot_pT(220, 325, "R1234yf_pT.png")
R1234ze = PR74(Tc=382.45, pc=3.6349, omega=0.313)
R1234ze.plot_pT(220, 325, "R1234ze(E)_pT.png")

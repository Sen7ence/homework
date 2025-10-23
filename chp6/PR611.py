import numpy as np
import matplotlib.pyplot as plt
import os

# 使用 Times New Roman 作为 matplotlib 全局字体
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14  # 增大全局字体
plt.rcParams["axes.linewidth"] = 1.5  # 增粗坐标轴


class PR611:
    def __init__(self, Tc1, pc1, omega1, Tc2, pc2, omega2, x1, kij):
        self.Tc1 = Tc1
        self.pc1 = pc1 * 1e6
        self.omega1 = omega1
        self.Tc2 = Tc2
        self.pc2 = pc2 * 1e6
        self.omega2 = omega2
        self.x1 = x1
        self.x2 = 1 - x1
        self.kij = kij

    R = 8.314462618  # J/(mol·K)

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
        return a1, a2, a, b1, b2, b, da

    # 计算A和B
    def AB(self, T, p):
        a1, a2, a, b1, b2, b, da = self.params(T)
        A = a * p * 1e6 / (self.R * T) ** 2
        B = b * p * 1e6 / (self.R * T)
        return A, B

    # 计算C2，C1，C0
    def C(self, T, p):
        A, B = self.AB(T, p)
        C2 = -(1 - B)
        C1 = A - 3 * B**2 - 2 * B
        C0 = -(A * B - B**2 - B**3)
        return C2, C1, C0

    # 计算压缩因子Z
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
    # 气相
    def vg(self, T, p):
        Zg = self.Zg(T, p)
        vg = Zg * self.R * T / (p * 1e6)
        return vg  # m³/mol

    # 计算逸度系数
    # 气相
    def phi_g(self, T, p):
        Zg = self.Zg(T, p)
        a1, a2, a, b1, b2, b, da = self.params(T)
        A, B = self.AB(T, p)
        phi_g1 = np.exp(
            (b1 / b) * (Zg - 1)
            - np.log(Zg - B)
            - A
            / (B * np.sqrt(8))
            * (
                2 * (self.x2 * (1 - self.kij) * (a1 * a2) ** 0.5 + self.x1 * a1) / a
                - b1 / b
            )
            * np.log((Zg + 2.414 * B) / (Zg - 0.414 * B))
        )
        phi_g2 = np.exp(
            (b2 / b) * (Zg - 1)
            - np.log(Zg - B)
            - A
            / (B * np.sqrt(8))
            * (
                2 * (self.x1 * (1 - self.kij) * (a1 * a2) ** 0.5 + self.x2 * a2) / a
                - b2 / b
            )
            * np.log((Zg + 2.414 * B) / (Zg - 0.414 * B))
        )
        return phi_g1, phi_g2

    # 计算逸度
    # 气相
    def f_g(self, T, p):  # MPa
        phi_g1, phi_g2 = self.phi_g(T, p)
        f_g1 = self.x1 * phi_g1 * p
        f_g2 = self.x2 * phi_g2 * p
        return f_g1, f_g2

    # 绘制溶液气相f-T、phi-T图
    def plot_fT(self, fluid_name, p, T_min, T_max, nT=11):
        T_grid = np.linspace(T_min, T_max, nT)  # 温度网格
        phi_grid1 = np.empty_like(T_grid)  # 组分1逸度系数网格
        phi_grid2 = np.empty_like(T_grid)  # 组分2逸度系数网格
        f_grid1 = np.empty_like(T_grid)  # 组分1逸度网格
        f_grid2 = np.empty_like(T_grid)  # 组分2逸度网格

        base_dir = os.path.dirname(os.path.abspath(__file__))
        fig_dir = os.path.join(base_dir, "figs")
        os.makedirs(fig_dir, exist_ok=True)

        # 计算逸度系数和逸度
        for i, T in enumerate(T_grid):
            phi_g1, phi_g2 = self.phi_g(T, p)
            f_g1, f_g2 = self.f_g(T, p)
            phi_grid1[i] = phi_g1
            phi_grid2[i] = phi_g2
            f_grid1[i] = f_g1
            f_grid2[i] = f_g2

        # T-phi
        fig_phi = plt.figure(figsize=(8, 6))
        ax_phi = fig_phi.add_subplot(1, 1, 1)
        ax_phi.plot(T_grid, phi_grid1, "b-o", linewidth=2, label="R290", markersize=6)
        ax_phi.plot(T_grid, phi_grid2, "r-s", linewidth=2, label="R600a", markersize=6)
        ax_phi.set_xlabel(r"$T$ (K)", fontsize=12)
        ax_phi.set_ylabel(r"$\hat{\phi}$", fontsize=12)
        # ax_phi.set_title(
        #    rf"{fluid_name} $\hat{{\phi}}$ at $p$ = {p:.1f} MPa", fontsize=12
        # )
        ax_phi.grid(True, linestyle="--", alpha=0.7)
        ax_phi.legend(loc="best", frameon=True, fancybox=True, framealpha=0.9)
        fig_phi.tight_layout()
        savepath_phi = os.path.join(fig_dir, f"{fluid_name}_phi.png")
        fig_phi.savefig(savepath_phi, dpi=600, bbox_inches="tight", transparent=False)
        plt.close(fig_phi)

        # T-f
        fig_f = plt.figure(figsize=(8, 6))
        ax_f = fig_f.add_subplot(111)
        ax_f.plot(T_grid, f_grid1, "b-o", linewidth=2, label="R290", markersize=6)
        ax_f.plot(T_grid, f_grid2, "r-s", linewidth=2, label="R600a", markersize=6)
        ax_f.set_xlabel(r"$T$ (K)", fontsize=12)
        ax_f.set_ylabel(r"$\hat{f}$ (MPa)", fontsize=12)
        # ax_f.set_title(rf"{fluid_name} $\hat{{f}}$ at $p$ = {p:.1f} MPa", fontsize=12)
        ax_f.grid(True, linestyle="--", alpha=0.7)
        ax_f.legend(loc="best", frameon=True, fancybox=True, framealpha=0.9)
        fig_f.tight_layout()
        savepath_f = os.path.join(fig_dir, f"{fluid_name}_f.png")
        fig_f.savefig(savepath_f, dpi=600, bbox_inches="tight", transparent=False)
        plt.close(fig_f)


R290R600a = PR611(
    Tc1=369.89,
    pc1=4.2512,
    omega1=0.1521,
    Tc2=407.81,
    pc2=3.629,
    omega2=0.184,
    x1=0.5,
    kij=0.064,
)

R290R600a.plot_fT("R290R600a", 1.4, 350, 450, 11)

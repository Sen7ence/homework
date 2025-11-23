import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14
plt.rcParams["axes.linewidth"] = 1.5


class PR75:
    def __init__(self, Tc1, pc1, omega1, Tc2, pc2, omega2, kij):
        self.Tc1 = float(Tc1)
        self.pc1 = float(pc1) * 1e6
        self.omega1 = float(omega1)
        self.Tc2 = float(Tc2)
        self.pc2 = float(pc2) * 1e6
        self.omega2 = float(omega2)
        self.kij = float(kij)

    R = 8.314462618  # J/(mol·K)

    # 计算 a_i, b_i 以及混合 a, b；mu1 为该相中组分1摩尔分数
    def params(self, T, mu1):
        kappa1 = 0.37464 + 1.54226 * self.omega1 - 0.26992 * self.omega1**2
        kappa2 = 0.37464 + 1.54226 * self.omega2 - 0.26992 * self.omega2**2
        Tr1 = T / self.Tc1
        Tr2 = T / self.Tc2
        alpha1 = (1 + kappa1 * (1 - np.sqrt(Tr1))) ** 2
        alpha2 = (1 + kappa2 * (1 - np.sqrt(Tr2))) ** 2
        a1 = 0.45724 * (self.R * self.Tc1) ** 2 / self.pc1 * alpha1
        a2 = 0.45724 * (self.R * self.Tc2) ** 2 / self.pc2 * alpha2
        b1 = 0.07780 * self.R * self.Tc1 / self.pc1
        b2 = 0.07780 * self.R * self.Tc2 / self.pc2
        a = (
            (mu1**2) * a1
            + (1 - mu1) ** 2 * a2
            + 2 * mu1 * (1 - mu1) * np.sqrt(a1 * a2) * (1 - self.kij)
        )
        b = mu1 * b1 + (1 - mu1) * b2
        return a1, a2, a, b1, b2, b

    # 计算 A 与 B
    def AB(self, T, p, mu1):
        a1, a2, a, b1, b2, b = self.params(T, mu1)
        p_Pa = p * 1e6
        A = a * p_Pa / (self.R**2 * T**2)
        B = b * p_Pa / (self.R * T)
        return A, B

    # PR 三次的系数
    def C(self, T, p, mu1):
        A, B = self.AB(T, p, mu1)
        C2 = B - 1.0
        C1 = A - 3.0 * B**2 - 2.0 * B
        C0 = -(A * B - B**2 - B**3)
        return C2, C1, C0

    # 液相 Z
    def Zl(self, T, p, mu1):
        C2, C1, C0 = self.C(T, p, mu1)
        Z = 1.0e-3
        for _ in range(10000):
            f = Z**3 + C2 * Z**2 + C1 * Z + C0
            df = 3.0 * Z**2 + 2.0 * C2 * Z + C1
            Z_new = Z - f / df
            if abs(Z_new - Z) < 1e-6:
                Z = Z_new
                break
            Z = Z_new
        return Z

    # 气相 Z
    def Zg(self, T, p, mu1):
        C2, C1, C0 = self.C(T, p, mu1)
        Z = 1.1
        for _ in range(10000):
            f = Z**3 + C2 * Z**2 + C1 * Z + C0
            df = 3.0 * Z**2 + 2.0 * C2 * Z + C1
            Z_new = Z - f / df
            if abs(Z_new - Z) < 1e-6:
                Z = Z_new
                break
            Z = Z_new
        return Z

    # 气相逸度系数（相组成取 mu1=y1）
    def phi_g(self, T, p, mu1):
        Zg = self.Zg(T, p, mu1)
        a1, a2, a, b1, b2, b = self.params(T, mu1)
        A, B = self.AB(T, p, mu1)
        # 组分1
        lnphi1 = (
            b1 / b * (Zg - 1.0)
            - np.log(Zg - B)
            - A
            / (2.0 * np.sqrt(2.0) * B)
            * (
                2.0 * ((1 - mu1) * (1 - self.kij) * np.sqrt(a1 * a2) + mu1 * a1) / a
                - b1 / b
            )
            * np.log((Zg + 2.414 * B) / (Zg - 0.414 * B))
        )
        # 组分2
        lnphi2 = (
            b2 / b * (Zg - 1.0)
            - np.log(Zg - B)
            - A
            / (2.0 * np.sqrt(2.0) * B)
            * (
                2.0 * (mu1 * (1 - self.kij) * np.sqrt(a1 * a2) + (1 - mu1) * a2) / a
                - b2 / b
            )
            * np.log((Zg + 2.414 * B) / (Zg - 0.414 * B))
        )
        return np.exp(lnphi1), np.exp(lnphi2)

    # 液相逸度系数（相组成取 mu1=x1）
    def phi_l(self, T, p, mu1):
        Zl = self.Zl(T, p, mu1)
        a1, a2, a, b1, b2, b = self.params(T, mu1)
        A, B = self.AB(T, p, mu1)

        lnphi1 = (
            b1 / b * (Zl - 1.0)
            - np.log(Zl - B)
            - A
            / (2.0 * np.sqrt(2.0) * B)
            * (
                2.0 * ((1 - mu1) * (1 - self.kij) * np.sqrt(a1 * a2) + mu1 * a1) / a
                - b1 / b
            )
            * np.log((Zl + 2.414 * B) / (Zl - 0.414 * B))
        )

        lnphi2 = (
            b2 / b * (Zl - 1.0)
            - np.log(Zl - B)
            - A
            / (2.0 * np.sqrt(2.0) * B)
            * (
                2.0 * (mu1 * (1 - self.kij) * np.sqrt(a1 * a2) + (1 - mu1) * a2) / a
                - b2 / b
            )
            * np.log((Zl + 2.414 * B) / (Zl - 0.414 * B))
        )
        return np.exp(lnphi1), np.exp(lnphi2)

    def plot_Tx(self, p, T0):
        T_list = []
        x_bub_list = []  # 泡点：x1 vs T
        y_dew_list = []  # 露点：y1 vs T

        # y1 扫描
        y1_values = np.linspace(0.0, 1.0, 1001)

        for y1 in y1_values:
            y1 = float(y1)
            y2 = 1.0 - y1

            # 初值：液相 x 猜 0.1/0.9；T 从 T0 逐步增加
            T = float(T0)
            x1 = 0.1
            s = 0.0  # s = Σ k_i y_i

            # 外层调温：使 s → 1（露点判据）
            # 注：气相侧 φ^v 用 y1；液相侧 φ^l 用 x1
            it_guard = 0
            while abs(s - 1.0) >= 1e-3:
                T += 0.1
                phi_g1, phi_g2 = self.phi_g(T, p, y1)  # 气相
                phi_l1, phi_l2 = self.phi_l(T, p, x1)  # 液相

                k1 = phi_g1 / phi_l1  # = 1/K1
                k2 = phi_g2 / phi_l2  # = 1/K2

                denom = k1 * y1 + k2 * y2
                if denom <= 1e-16:
                    break

                x1 = (k1 * y1) / denom
                s_prev = s
                s = denom

                # 细化循环：仅重算液相侧
                inner_guard = 0
                while abs(s - s_prev) > 1e-6:
                    phi_l1, phi_l2 = self.phi_l(T, p, x1)
                    k1 = phi_g1 / phi_l1
                    k2 = phi_g2 / phi_l2
                    denom = k1 * y1 + k2 * y2
                    if denom <= 1e-16:
                        break
                    x1 = (k1 * y1) / denom
                    s_prev = s
                    s = denom

                    inner_guard += 1
                    if inner_guard > 2000:  # 防止极端情况
                        break

                it_guard += 1
                if it_guard > 20000:  # 防止极端情况
                    break

            # 收集边界点
            T_list.append(T)
            x_bub_list.append(x1)  # 泡点边界（液相组成）
            y_dew_list.append(y1)  # 露点边界（气相组成）

        fig, ax = plt.subplots(1, figsize=(8, 6))
        ax.scatter(x_bub_list, T_list, s=3, label=r"bubble: $T$–$x$")
        ax.scatter(y_dew_list, T_list, s=3, label=r"dew: $T$–$y$")
        ax.set_xlim(0.0, 1.0)
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.grid(True)

        ax.set_xlabel(r"$x_\mathrm{R600a}$")
        ax.set_ylabel(r"$T$ (K)")
        ax.set_ylim(225, 265)

        base_dir = os.path.dirname(os.path.abspath(__file__))
        fig_dir = os.path.join(base_dir, "figs")
        os.makedirs(fig_dir, exist_ok=True)
        filename = f"{p:.3f}MPa.png"
        savepath = os.path.join(fig_dir, filename)
        fig.savefig(savepath, dpi=600, bbox_inches="tight", transparent=False)
        plt.close(fig)


R290R600a = PR75(
    Tc1=369.89,
    pc1=4.2512,
    omega1=0.1521,  # R290
    Tc2=407.81,
    pc2=3.629,
    omega2=0.184,  # R600a
    kij=0.01,
)
R290R600a.plot_Tx(p=0.1, T0=215.0)
# R290R600a.plot_Tx(p=1.0, T0=290.0)

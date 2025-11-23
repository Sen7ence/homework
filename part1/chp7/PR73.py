import numpy as np
import matplotlib.pyplot as plt
import os

# 使用 Times New Roman 作为 matplotlib 全局字体
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14  # 增大全局字体
plt.rcParams["axes.linewidth"] = 1.5  # 增粗坐标轴


class PR73:
    def __init__(self, Tc, pc, A1, A2, A3):
        self.Tc = Tc  # K
        self.pc = pc * 1e6  # Pa，输入MPa
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3

    def psat(self, T):
        Tr = T / self.Tc
        p = self.pc * (Tr) ** (
            self.A1 + self.A2 * (1 - Tr) ** 1.89 + self.A3 * (1 - Tr) ** 5.67
        )
        return p / 1e6  # 输出MPa

    def plot_psat(self, fluid_name, Tmin, Tmax, nt):
        T_grid = np.linspace(Tmin, Tmax, nt)
        psat_grid = np.zeros(nt)
        for i, T in enumerate(T_grid):
            psat_grid[i] = self.psat(T)

        base_dir = os.path.dirname(os.path.abspath(__file__))
        fig_dir = os.path.join(base_dir, "figs")
        os.makedirs(fig_dir, exist_ok=True)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(T_grid, psat_grid, "b-", linewidth=2)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$T$(K)", fontsize=14)
        ax.set_ylabel(r"$p_\mathrm{sat}$(MPa)", fontsize=14)
        ax.grid(True)

        # 保存图像
        savepath = os.path.join(fig_dir, f"{fluid_name}.png")
        fig.savefig(savepath, dpi=600, bbox_inches="tight", transparent=False)
        plt.close(fig)


CH4 = PR73(190.551, 4.5992, 5.87304544, 6.23280143, 13.0721578)
CH4.plot_psat("甲烷", 80, 180, 100)
C2H6 = PR73(305.33, 4.8717, 6.30717658, 7.47042131, 17.0958137)
C2H6.plot_psat("乙烷", 200, 300, 100)
C3H8 = PR73(369.80, 4.239, 6.50580501, 8.6776247, 18.0116214)
C3H8.plot_psat("丙烷", 200, 350, 150)
C4H10 = PR73(425.2, 3.8, 6.81692028, 8.77671813, 23.7680492)
C4H10.plot_psat("丁烷", 200, 350, 150)

import numpy as np
import matplotlib.pyplot as plt
import os

# 使用 Times New Roman 作为 matplotlib 全局字体
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14  # 增大全局字体
plt.rcParams["axes.linewidth"] = 1.5  # 增粗坐标轴


class PR:
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
        return a1, a2, a, b1, b2, b, da

    # 计算A和B
    def AB(self, T, p):
        a1, a2, a, b1, b2, b, da = self.params(T)
        A = a * p * 1e6 / (self.R * T) ** 2
        B = b * p * 1e6 / (self.R * T)
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
        vl = Zl * self.R * T / (p * 1e6)
        return vl  # m³/mol

    # 气相
    def vg(self, T, p):
        Zg = self.Zg(T, p)
        vg = Zg * self.R * T / (p * 1e6)
        return vg  # m³/mol

    def plot_Tv(
        self,
        fluid_name,  # 流体名称
        p,  # 压力 MPa
        Tsat,  # 饱和温度 K
        T_min,  # 温度范围最小值 K
        T_max,  # 温度范围最大值 K
        nT=220,  # 温度点数
    ):
        T_grid = np.linspace(T_min, T_max, nT)  # 温度网格
        v_grid = np.empty_like(T_grid)  # 比体积网格
        # 计算比体积
        for i, T in enumerate(T_grid):
            if T < Tsat:
                v_grid[i] = self.vl(T, p) / (self.x1 * self.M1 + self.x2 * self.M2)
            elif T > Tsat:
                v_grid[i] = self.vg(T, p) / (self.x1 * self.M1 + self.x2 * self.M2)
            else:
                v_grid[i] = (
                    0.5
                    * (self.vl(T, p) + self.vg(T, p))
                    / (self.x1 * self.M1 + self.x2 * self.M2)
                )
        fig, ax = plt.subplots(1, figsize=(8, 6))  # 创建图像和坐标轴
        # 主曲线
        ax.plot(
            v_grid,
            T_grid,
            linewidth=2,
            label=fluid_name,
        )
        xmin, xmax = np.nanmin(v_grid), np.nanmax(v_grid)
        # Tsat 虚线
        ax.hlines(Tsat, xmin, xmax, linestyles="--", label=r"$T_{\mathrm{sat}}$")
        # 标注 Tsat
        yt = list(ax.get_yticks())
        # 加入Tsat并排序
        if not any(abs(t - Tsat) < 1e-8 for t in yt):
            yt.append(Tsat)
        yt = np.array(sorted(yt))
        # 生成刻度标签：对 Tsat 使用仅数值标签（两位小数），其它刻度保留数字格式（根据范围选择小数位）
        deltaT = T_grid.max() - T_grid.min()
        labels = []
        for t in yt:
            if abs(t - Tsat) < 1e-8 or abs(t - Tsat) < 1e-6 * max(1.0, deltaT):
                labels.append(f"{Tsat:.2f}")
            else:
                # 根据温度范围决定格式，避免过多小数
                if deltaT > 50:
                    labels.append(f"{t:.0f}")
                else:
                    labels.append(f"{t:.2f}")
        ax.set_yticks(yt)
        ax.set_yticklabels(labels)
        # 轴标签
        ax.set_xlabel(r"$v$ (m³/kg)")
        ax.set_ylabel(r"$T$ (K)")
        # 标题
        # ax.set_title(f"{fluid_name}  $T$–$v$ at $p$ = {p:.1f} MPa")
        ax.grid(True)
        ax.set_xscale("log")  # 使用对数刻度
        ax.legend(loc="upper left", frameon=True, fancybox=True, framealpha=0.9)

        # 固定保存路径为脚本同目录下的 figs 文件夹
        base_dir = os.path.dirname(os.path.abspath(__file__))
        fig_dir = os.path.join(base_dir, "figs")
        os.makedirs(fig_dir, exist_ok=True)

        # 文件名固定为"流体名称.png"
        filename = f"{fluid_name}.png"
        savepath = os.path.join(fig_dir, filename)

        # 保存图像，固定参数
        fig.savefig(savepath, dpi=600, bbox_inches="tight", transparent=False)
        plt.close(fig)

    # 计算焓的余函数
    # 液相
    def h_res_l(self, T, p):
        a1, a2, a, b1, b2, b, da = self.params(T)
        Zl = self.Zl(T, p)
        vl = self.vl(T, p)
        hr_l = (T * da - a) / (b * np.sqrt(8)) * np.log(
            (vl - 0.414 * b) / (vl + 2.414 * b)
        ) + self.R * T * (1 - Zl)
        return hr_l

    # 气相
    def h_res_g(self, T, p):
        a1, a2, a, b1, b2, b, da = self.params(T)
        Zg = self.Zg(T, p)
        vg = self.vg(T, p)
        hr_g = (T * da - a) / (b * np.sqrt(8)) * np.log(
            (vg - 0.414 * b) / (vg + 2.414 * b)
        ) + self.R * T * (1 - Zg)
        return hr_g

    # 计算熵的余函数
    # 液相
    def s_res_l(self, T, p):
        a1, a2, a, b1, b2, b, da = self.params(T)
        vl = self.vl(T, p)
        sr_l = (
            -self.R * np.log((vl - b) / vl)
            - self.R * np.log(vl / (self.R * T / (p * 1e6)))
            + da / (b * np.sqrt(8)) * np.log((vl - 0.414 * b) / (vl + 2.414 * b))
        )
        return sr_l

    # 气相
    def s_res_g(self, T, p):
        a1, a2, a, b1, b2, b, da = self.params(T)
        vg = self.vg(T, p)
        sr_g = (
            -self.R * np.log((vg - b) / vg)
            - self.R * np.log(vg / (self.R * T / (p * 1e6)))
            + da / (b * np.sqrt(8)) * np.log((vg - 0.414 * b) / (vg + 2.414 * b))
        )
        return sr_g

    # 计算c_p积分
    def cp(self, T, A, B, C, D):
        cp = (
            A * (T - 273.15)
            + B / 2 * (T**2 - 273.15**2)
            + C / 3 * (T**3 - 273.15**3)
            + D / 4 * (T**4 - 273.15**4)
        )
        return cp

    # 计算c_p/T积分
    def cpT(self, T, A, B, C, D):
        cp = (
            A * np.log(T / 273.15)
            + B * (T - 273.15)
            + C / 2 * (T**2 - 273.15**2)
            + D / 3 * (T**3 - 273.15**3)
        )
        return cp

    # 计算焓和熵
    # 液相
    def h_l(self, T, A, B, C, D, p):
        h_r_ps_0 = self.h_res_l(273.15, self.ps0)
        cp0 = self.cp(T, A, B, C, D)
        h_res_l = self.h_res_l(T, p)
        hl = (
            200 * 1e3
            + cp0
            + (h_r_ps_0 - h_res_l) / (self.x1 * self.M1 + self.x2 * self.M2)
        )  # J/kg
        return hl

    def s_l(self, T, A, B, C, D, p):
        s_r_ps_0 = self.s_res_l(273.15, self.ps0)
        cpT = self.cpT(T, A, B, C, D)
        sr_l = self.s_res_l(T, p)
        sl = (
            1e3
            + cpT
            + (s_r_ps_0 - self.R * np.log(p / self.ps0) - sr_l)
            / (self.x1 * self.M1 + self.x2 * self.M2)
        )  # J/(kg*K)
        return sl

    # 气相
    def h_g(self, T, A, B, C, D, p):
        h_r_ps_0 = self.h_res_l(273.15, self.ps0)  # 使用液相作为基准
        cp0 = self.cp(T, A, B, C, D)
        h_res_g = self.h_res_g(T, p)
        hg = (
            200 * 1e3
            + cp0
            + (h_r_ps_0 - h_res_g) / (self.x1 * self.M1 + self.x2 * self.M2)
        )  # J/kg
        return hg

    def s_g(self, T, A, B, C, D, p):
        s_r_ps_0 = self.s_res_l(273.15, self.ps0)  # 使用液相作为基准
        cpT = self.cpT(T, A, B, C, D)
        sr_g = self.s_res_g(T, p)
        sg = (
            1e3
            + cpT
            + (s_r_ps_0 - self.R * np.log(p / self.ps0) - sr_g)
            / (self.x1 * self.M1 + self.x2 * self.M2)
        )  # J/(kg*K)
        return sg

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

    # 液相
    def phi_l(self, T, p):
        Zl = self.Zl(T, p)
        a1, a2, a, b1, b2, b, da = self.params(T)
        A, B = self.AB(T, p)
        phi_l1 = np.exp(
            (b1 / b) * (Zl - 1)
            - np.log(Zl - B)
            - A
            / (B * np.sqrt(8))
            * (
                2 * (self.x2 * (1 - self.kij) * (a1 * a2) ** 0.5 + self.x1 * a1) / a
                - b1 / b
            )
            * np.log((Zl + 2.414 * B) / (Zl - 0.414 * B))
        )
        phi_l2 = np.exp(
            (b2 / b) * (Zl - 1)
            - np.log(Zl - B)
            - A
            / (B * np.sqrt(8))
            * (
                2 * (self.x1 * (1 - self.kij) * (a1 * a2) ** 0.5 + self.x2 * a2) / a
                - b2 / b
            )
            * np.log((Zl + 2.414 * B) / (Zl - 0.414 * B))
        )
        return phi_l1, phi_l2

    # 计算逸度
    # 气相
    def f_g(self, T, p):  # MPa
        phi_g1, phi_g2 = self.phi_g(T, p)
        f_g1 = self.x1 * phi_g1 * p
        f_g2 = self.x2 * phi_g2 * p
        return f_g1, f_g2

    # 液相
    def f_l(self, T, p):  # MPa
        phi_l1, phi_l2 = self.phi_l(T, p)
        f_l1 = self.x1 * phi_l1 * p
        f_l2 = self.x2 * phi_l2 * p
        return f_l1, f_l2

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

    # p-T相图
    def plot_pT(self, fluid_name, T_min, T_max, nT=220):
        T_grid = np.linspace(T_min, T_max, nT)  # 温度网格
        p_grid = np.empty_like(T_grid)  # 压力网格
        # 计算饱和压力
        for i, T in enumerate(T_grid):
            # 使用二分法求解饱和压力
            p_low = 0.001  # MPa
            p_high = 10.0  # MPa
            for _ in range(100):
                p_mid = (p_low + p_high) / 2
                phi_l1, phi_l2 = self.phi_l(T, p_mid)
                phi_g1, phi_g2 = self.phi_g(T, p_mid)
                f_l1 = self.x1 * phi_l1 * p_mid
                f_l2 = self.x2 * phi_l2 * p_mid
                f_g1 = self.x1 * phi_g1 * p_mid
                f_g2 = self.x2 * phi_g2 * p_mid
                if (f_l1 - f_g1) > 0:
                    p_high = p_mid
                else:
                    p_low = p_mid
                if abs(p_high - p_low) < 1e-6:
                    break
            p_grid[i] = p_mid
        fig, ax = plt.subplots(1, figsize=(8, 6))  # 创建图像和坐标轴
        # 主曲线
        ax.plot(
            T_grid,
            p_grid,
            linewidth=2,
            label=fluid_name,
        )
        # 轴标签
        ax.set_xlabel(r"$T$ (K)")
        ax.set_ylabel(r"$p$ (MPa)")
        # 标题
        # ax.set_title(f"{fluid_name}  $p$–$T$ Saturation Curve")
        ax.grid(True)
        ax.legend(loc="upper left", frameon=True, fancybox=True, framealpha=0.9)

        # 固定保存路径为脚本同目录下的 figs 文件夹
        base_dir = os.path.dirname(os.path.abspath(__file__))
        fig_dir = os.path.join(base_dir, "figs")
        os.makedirs(fig_dir, exist_ok=True)

        # 文件名固定为"流体名称_pT.png"
        filename = f"{fluid_name}_pT.png"
        savepath = os.path.join(fig_dir, filename)
        # 保存图像，固定参数
        fig.savefig(savepath, dpi=600, bbox_inches="tight", transparent=False)
        plt.close(fig)


if __name__ == "__main__":
    # 3-10
    # R290 = PR(
    #     Tc1=369.89,  # K
    #     pc1=4.2512,  # MPa
    #     omega1=0.1521,  # 无量纲
    #     M1=44.096,  # g/mol
    #     x1=1.0,
    #     Tc2=369.89,
    #     pc2=4.2512,
    #     omega2=0.1521,
    #     M2=44.096,
    #     kij=0.064,
    #     ps0=0.47446,
    # )

    # R600a = PR(
    #     Tc1=407.81,
    #     pc1=3.629,
    #     omega1=0.184,
    #     M1=58.122,  # R600a
    #     x1=1.0,
    #     Tc2=407.81,
    #     pc2=3.629,
    #     omega2=0.184,
    #     M2=58.122,  # R600a
    #     kij=0.0,
    #     ps0=0.15696,
    # )

    # R290.plot_Tv("R290", 1.4, 317.86, 200, 450)
    # R600a.plot_Tv("R600a", 0.6, 314.12, 200, 450)

    # 3-13
    # R134a = PR(
    #     Tc1=374.21,
    #     pc1=4.0593,
    #     omega1=0.326,
    #     M1=102.03,  # R134a
    #     x1=1.0,
    #     Tc2=374.21,
    #     pc2=4.0593,
    #     omega2=0.326,
    #     M2=102.03,  # R134a
    #     kij=0.0,
    #     ps0=0.57245,
    # )
    # R1234yf = PR(
    #     Tc1=367.85,
    #     pc1=3.3822,
    #     omega1=0.276,
    #     M1=114.04,  # R1234yf
    #     x1=1.0,
    #     Tc2=367.85,
    #     pc2=3.3822,
    #     omega2=0.276,
    #     M2=114.04,  # R1234yf
    #     kij=0.0,
    #     ps0=0.42483,
    # )
    # R1234ze = PR(
    #     Tc1=382.75,
    #     pc1=3.6349,
    #     omega1=0.313,
    #     M1=114.04,  # R1234ze
    #     x1=1.0,
    #     Tc2=382.75,
    #     pc2=3.6349,
    #     omega2=0.313,
    #     M2=114.04,  # R1234ze
    #     kij=0.0,
    #     ps0=0.49314,
    # )

    # print(R134a.vg(308.15, 0.1) / (R134a.x1 * R134a.M1 + R134a.x2 * R134a.M2))
    # print(R1234yf.vg(308.15, 0.1) / (R1234yf.x1 * R1234yf.M1 + R1234yf.x2 * R1234yf.M2))
    # print(R1234ze.vg(308.15, 0.1) / (R1234ze.x1 * R1234ze.M1 + R1234ze.x2 * R1234ze.M2))

    # 3-15
    # kij = 0.064, p=0.1 MPa情况下
    # R290R600a = PR(
    #     Tc1=369.89,
    #     pc1=4.2512,
    #     omega1=0.1521,
    #     M1=44.096,  # R290
    #     x1=0.5,
    #     Tc2=407.81,
    #     pc2=3.629,
    #     omega2=0.184,
    #     M2=58.122,  # R600a
    #     kij=0.064,
    #     ps0=0.32979,
    # )

    # print(
    #     R290R600a.vg(300, 0.1)
    #     / (R290R600a.x1 * R290R600a.M1 + R290R600a.x2 * R290R600a.M2)
    # )

    # 4-13
    # R290 = PR(
    #     Tc1=369.89,  # K
    #     pc1=4.2512,  # MPa
    #     omega1=0.1521,  # 无量纲
    #     M1=44.096,  # g/mol
    #     x1=1.0,
    #     Tc2=369.89,
    #     pc2=4.2512,
    #     omega2=0.1521,
    #     M2=44.096,
    #     kij=0.064,
    #     ps0=0.47446,
    # )

    # R600a = PR(
    #     Tc1=407.81,
    #     pc1=3.629,
    #     omega1=0.184,
    #     M1=58.122,  # R600a
    #     x1=1.0,
    #     Tc2=407.81,
    #     pc2=3.629,
    #     omega2=0.184,
    #     M2=58.122,  # R600a
    #     kij=0.0,
    #     ps0=0.15696,
    # )
    # # 300K下结果
    # print(R290.h_l(300, -95.80, 6.945, -3.597 * 1e-3, 7.290 * 1e-7, 1.4))
    # print(R290.s_l(300, -95.80, 6.945, -3.597 * 1e-3, 7.290 * 1e-7, 1.4))
    # print(R600a.h_l(300, -23.91, 6.605, -3.176 * 1e-3, 4.981 * 1e-7, 0.6))
    # print(R600a.s_l(300, -23.91, 6.605, -3.176 * 1e-3, 4.981 * 1e-7, 0.6))

    # 4-15
    # R290R600a = PR(
    #     Tc1=369.89,
    #     pc1=4.2512,
    #     omega1=0.1521,
    #     M1=44.096,  # R290
    #     x1=0.5,
    #     Tc2=407.81,
    #     pc2=3.629,
    #     omega2=0.184,
    #     M2=58.122,  # R600a
    #     kij=0.064,
    #     ps0=0.32979,
    # )
    # # 300K下结果
    # print(R290R600a.h_l(260, -59.81, 6.775, -3.386 * 1e-3, 6.135 * 1e-7, 1))
    # print(R290R600a.s_l(260, -59.81, 6.775, -3.386 * 1e-3, 6.135 * 1e-7, 1))

    # 6-11
    # R290R600a = PR(
    #     Tc1=369.89,
    #     pc1=4.2512,
    #     omega1=0.1521,
    #     M1=44.096,  # R290
    #     x1=0.5,
    #     Tc2=407.81,
    #     pc2=3.629,
    #     omega2=0.184,
    #     M2=58.122,  # R600a
    #     kij=0.064,
    #     ps0=0.32979,
    # )
    # R290R600a.plot_fT("R290R600a", 1.4, 350, 450, 11)

    # 7-4
    R290 = PR(
        Tc1=369.89,  # K
        pc1=4.2512,  # MPa
        omega1=0.1521,  # 无量纲
        M1=44.096,  # g/mol
        x1=1.0,
        Tc2=369.89,
        pc2=4.2512,
        omega2=0.1521,
        M2=44.096,
        kij=0.064,
        ps0=0.47446,
    )
    R290.plot_pT("R290", 200, 450, 220)
    # 7-5

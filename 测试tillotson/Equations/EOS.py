"""
    状态方程
    状态方程负责计算粒子压力、声速
    有时还需要计算温度
"""

from pysph.sph.equation import Equation
from math import exp
import config


class TaitEOS(Equation):
    """
        最简单的EOS
    """
    def __init__(self, dest, sources, rho0, c0, gamma, p0=0.0):
        self.rho0 = rho0
        self.rho01 = 1.0 / rho0
        self.c0 = c0
        self.gamma = gamma
        self.gamma1 = 0.5 * (gamma - 1.0)
        self.B = rho0 * c0 * c0 / gamma
        self.p0 = p0

        super(TaitEOS, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho, d_p, d_cs):
        ratio = d_rho[d_idx] * self.rho01
        tmp = pow(ratio, self.gamma)

        d_p[d_idx] = self.p0 + self.B * (tmp - 1.0)
        d_cs[d_idx] = self.c0 * pow(ratio, self.gamma1)


class crashEOS(Equation):
    """
        一个固体的EOS，但不适用于高速撞击
    """
    def __init__(self, dest, sources, rho0, k1=45.4, k2=-138, k3=290, gamma=7, c0=2000):
        self.rho0 = rho0
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.rho01 = 1.0 / rho0
        self.gamma1 = 0.5 * (gamma - 1.0)
        self.c0 = c0

        super(crashEOS, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho, d_p, d_cs, d_T):
        miu = d_rho[d_idx] / self.rho0 - 1
        if miu > 0:
            d_p[d_idx] = miu * (self.k1 + miu * (self.k2 + miu * self.k3))
        else:
            d_p[d_idx] = miu * self.k1
        ratio = d_rho[d_idx] * self.rho01
        d_cs[d_idx] = self.c0 * pow(ratio, self.gamma1)
        d_T[d_idx] = 400


class TillotsonEOS(Equation):
    """
        最合适的EOS
        但数量级可能还需要根据实际情况调
    """

    def __init__(self, dest, sources):
        super(TillotsonEOS, self).__init__(dest, sources)

    def loop(self, d_p, d_idx, d_rho, d_e, d_cs, d_T):

        # 1erg/g = 1e4J/kg
        e = d_e[d_idx] * 1e4

        # 1g/cm2 = 1e3kg/m3
        rho = d_rho[d_idx] * 1e-3

        eta = rho / config.tillotson_rho0

        miu = eta - 1

        if e < config.tillotson_us or eta >= 1:
            d_p[d_idx] = (config.tillotson_a + config.tillotson_b / (
                    e / config.tillotson_e0 / eta ** 2 + 1)) * e * rho + config.tillotson_A * miu + config.tillotson_B * miu ** 2
        elif e > config.tillotson_us_:
            z = 1 / eta - 1
            ea = exp(-config.tillotson_alpha * z)
            eb = exp(-config.tillotson_beta * z ** 2)
            d_p[d_idx] = config.tillotson_a * e * rho + (config.tillotson_b * e * rho / (
                        e / config.tillotson_e0 / eta ** 2 + 1) + config.tillotson_B * miu * ea) * eb
        else:
            z = 1 / eta - 1
            ea = exp(-config.tillotson_alpha * z)
            eb = exp(-config.tillotson_beta * z ** 2)
            e1 = e - config.tillotson_us
            pe = config.tillotson_a * e1 * rho + (config.tillotson_b * e1 * rho / (
                        e1 / config.tillotson_e0 / eta ** 2 + 1) + config.tillotson_A * miu * ea) * eb
            e2 = config.tillotson_us_ - e
            pc = (config.tillotson_a + config.tillotson_b / (
                        e2 / config.tillotson_e0 / eta ** 2 + 1)) * e2 * rho + config.tillotson_A * miu + config.tillotson_B * miu ** 2
            d_p[d_idx] = (pe + pc) / (config.tillotson_us_ - config.tillotson_us)

        ratio = d_rho[d_idx] / config.tillotson_rho0
        d_cs[d_idx] = config.tillotson_c0 * pow(ratio, 7)
        # d_cs[d_idx] = 1e-4
        d_T[d_idx] = e / config.tillotson_cv

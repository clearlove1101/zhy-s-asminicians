"""
    状态方程
    状态方程负责计算粒子压力、声速
    有时还需要计算温度
"""

from pysph.sph.equation import Equation
from math import exp
import config


class TaitEOS(Equation):
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
    def __init__(self, dest, sources, rho0, k1=45.4, k2=-138, k3=290, gamma=7, c0=2000, cv=config.tillotson_cv):
        self.rho0 = rho0
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.rho01 = 1.0 / rho0
        self.gamma1 = 0.5 * (gamma - 1.0)
        self.c0 = c0
        self.cv = cv

        super(crashEOS, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho, d_p, d_cs, d_T, d_e):
        miu = d_rho[d_idx] / self.rho0 - 1
        if miu > 0:
            d_p[d_idx] = miu * (self.k1 + miu * (self.k2 + miu * self.k3))
            d_cs[d_idx] = ((self.k1 + 2 * miu * self.k2 + 3 * miu ** 2 * self.k3) * self.rho01) ** 0.5
        else:
            d_p[d_idx] = miu * self.k1
            d_cs[d_idx] = (self.k1 * self.rho01) ** 0.5
        # ratio = d_rho[d_idx] * self.rho01
        # d_cs[d_idx] = self.c0 * pow(ratio, self.gamma1)
        d_T[d_idx] = d_e[d_idx] / self.cv


class TillotsonEOS(Equation):

    def __init__(self,
                 dest, sources, rho0=config.tillotson_rho0, us=config.tillotson_us, a=config.tillotson_a,
                 b=config.tillotson_b, e0=config.tillotson_e0, A=config.tillotson_A, B=config.tillotson_B,
                 us_=config.tillotson_us_, alpha=config.tillotson_alpha, beta=config.tillotson_beta,
                 c0=config.tillotson_c0, cv=config.tillotson_cv
                 ):
        super(TillotsonEOS, self).__init__(dest, sources)
        self.rho0 = rho0
        self.us = us
        self.a = a
        self.b = b
        self.e0 = e0
        self.A = A
        self.B = B
        self.us_ = us_
        self.alpha = alpha
        self.beta = beta
        self.c0 = c0
        self.cv = cv


    def loop(self, d_p, d_idx, d_rho, d_e, d_cs, d_T):

        # 1erg/g = 1e4J/kg
        e = d_e[d_idx] * 1e4

        # 1g/cm3 = 1e3kg/m3
        rho = d_rho[d_idx] * 1e-3

        eta = rho / self.rho0

        miu = eta - 1

        if e < self.us or eta >= 1:

            w0 = e / self.e0 / eta ** 2 + 1.

            abedw = (self.a + self.b / w0) * e
            apbm = self.A + self.B * miu

            d_p[d_idx] = (abedw * rho + apbm * miu) * 0.1

            c2 = abedw + 2 * self.b * e ** 2 * rho / w0 ** 2 / self.rho0 / self.e0 / eta ** 3 + \
                 apbm / self.rho0

            d_cs[d_idx] = c2 ** 0.5 * 1e-2

        elif e > self.us_:
            z = 1 / eta - 1
            ea = exp(-self.alpha * z)
            betaz = self.beta * z
            eb = exp(-betaz * z)
            w0 = e / self.e0 / eta ** 2 + 1.
            berdw = self.b * e * rho / w0


            d_p[d_idx] = (self.a * e * rho + (berdw + self.A * miu * ea) * eb) * 0.1

            c2 = (self.a + self.b * eb / w0) * e +\
                2 * berdw * eb / self.rho0 / eta ** 2 * (e / self.e0 / w0 / eta - betaz) +\
                self.A * ea * eb / self.rho0 * (1. + self.alpha + 2 * betaz)

            d_cs[d_idx] = c2 ** 0.5 * 1e-2

        else:
            z = 1 / eta - 1
            w0 = e / self.e0 / eta ** 2 + 1.
            betaz = self.beta * z
            ea = exp(-self.alpha * z)
            eb = exp(-betaz * z)
            abedw = (self.a + self.b / w0) * e
            apbm = self.A + self.B * miu
            berdw = self.b * e * rho / w0

            pe = (self.a * e * rho + (berdw + self.A * miu * ea) * eb)
            pc = (abedw * rho + apbm * miu)

            d_p[d_idx] = ((pe * (e - self.us) + pc * (self.us_ - e)) / (self.us_ - self.us)) * 0.1

            c2_1 = abedw + 2 * self.b * e ** 2 * rho / w0 ** 2 / self.rho0 / self.e0 / eta ** 3 + \
                 apbm / self.rho0
            c2_2 = (self.a + self.b * eb / w0) * e +\
                2 * berdw * eb / self.rho0 / eta ** 2 * (e / self.e0 / w0 / eta - betaz) +\
                self.A * ea * eb / self.rho0 * (1. + self.alpha + 2 * betaz)
            c2 = ((c2_2 * (e - self.us) + c2_1 * (self.us_ - e)) / (self.us_ - self.us))
            d_cs[d_idx] = c2 ** 0.5 * 1e-2

        # d_cs[d_idx] = self.c0 * pow(eta, 3)
        # d_cs[d_idx] = 1e-4
        d_T[d_idx] = e / self.cv

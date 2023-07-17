from pysph.sph.equation import Equation
from math import exp
import config

class TillotsonEOS(Equation):
    """
        最合适的EOS
        但数量级可能还需要根据实际情况调
    """

    def __init__(self, dest, sources):
        super(TillotsonEOS, self).__init__(dest, sources)

    def loop(self, d_p, d_idx, d_rho, d_e, d_cs, d_T):

        # 1erg/g = 1e-4J/kg
        e = d_e[d_idx] * 1e4#SI换算成erg/g,所以为什么不用u，用e

        # 1g/cm3 = 1e3kg/m3
        rho = d_rho[d_idx] * 1e-3#SI换算成g/cm3


        eta = rho / config.tillotson_rho0

        miu = eta - 1

        if e < config.tillotson_us or eta >= 1:
            d_p[d_idx] = (config.tillotson_a + config.tillotson_b / (
                    e / config.tillotson_e0 / eta ** 2 + 1)) * e * rho
            + config.tillotson_A * miu + config.tillotson_B * miu ** 2
            d_p[d_idx] = d_p[d_idx] * 0.1#将达因/平方厘米换算成SI
        elif e > config.tillotson_us_:
            z = 1 / eta - 1
            ea = exp(-config.tillotson_alpha * z)
            eb = exp(-config.tillotson_beta * z ** 2)
            d_p[d_idx] = config.tillotson_a * e * rho + (config.tillotson_b * e * rho / (
                        e / config.tillotson_e0 / eta ** 2 + 1) + config.tillotson_B * miu * ea) * eb
            d_p[d_idx] = d_p[d_idx] * 0.1 
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
            d_p[d_idx] = d_p[d_idx] * 0.1
            
        d_cs[d_idx] = self.c0 * pow(eta, 3)
        # d_cs[d_idx] = 1e-4
        d_T[d_idx] = e / self.cv
'''
class TillotsonEOS(Equation):

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
            d_p[d_idx] = d_p[d_idx] * 0.1
        elif e > config.tillotson_us_:
            z = 1 / eta - 1
            ea = exp(-config.tillotson_alpha * z)
            eb = exp(-config.tillotson_beta * z ** 2)
            d_p[d_idx] = config.tillotson_a * e * rho + (config.tillotson_b * e * rho / (
                        e / config.tillotson_e0 / eta ** 2 + 1) + config.tillotson_B * miu * ea) * eb
            d_p[d_idx] = d_p[d_idx] * 0.1
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
            d_p[d_idx] = d_p[d_idx] * 0.1

        ratio = d_rho[d_idx] / config.tillotson_rho0
        d_cs[d_idx] = config.tillotson_c0 * pow(ratio, 7)
        # d_cs[d_idx] = 1e-4
        d_T[d_idx] = e / config.tillotson_cv
'''

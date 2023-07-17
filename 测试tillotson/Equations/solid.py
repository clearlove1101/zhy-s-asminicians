"""
    待实现ROCK等固体相关方程
"""
from pysph.sph.equation import Equation
from math import tanh
import config


class ROCK(Equation):
    def __init__(self, dest, sources):
        # kw = rock_parameters()
        self.Tm0 = config.rock_Tm0
        self.a = config.rock_a
        self.c = config.rock_c
        self.Xi = config.rock_Xi
        self.Epsilonfb = config.rock_Epsilonfb
        self.Pc = config.rock_Pc
        self.B = config.rock_B
        self.Yi0 = config.rock_Yi0
        self.miui = config.rock_miui
        self.Yim = config.rock_Yim
        self.Yd0 = config.rock_Yd0
        self.miud = config.rock_miud
        self.Ydm = config.rock_Ydm
        super(ROCK, self).__init__(dest, sources)

    def loop(self, d_idx, d_Y,
             d_T, d_p, d_dam, d_D
             ):
        dam = d_dam[d_idx]
        p = d_p[d_idx]
        Tm = self.Tm0 * (p / self.a + 1) ** (1 / self.c)
        YtDY = tanh(self.Xi * (Tm / d_T[d_idx] - 1))
        Epsilonf_ = self.B * (p - self.Pc)
        Epsilonf = Epsilonf_ if Epsilonf_ > self.Epsilonfb else self.Epsilonfb
        D_ = dam / Epsilonf
        if 0 < D_ < 1:
            D = D_
        elif D_ > 0:
            D = 1
        else:
            D = 0
        d_D[d_idx] = D
        tmp = self.miui * p
        Yi = self.Yi0 + 1 / (1 / tmp + 1 / (self.Yim - self.Yi0))
        Yd_ = self.Yd0 + self.miud * p
        Yd = Yd_ if Yd_ < self.Ydm else self.Ydm
        d_Y[d_idx] = (D * Yd + (1 - D) * Yi) * YtDY

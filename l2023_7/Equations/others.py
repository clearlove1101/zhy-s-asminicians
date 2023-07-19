"""
    其它方程
"""
from pysph.sph.equation import Equation


class InternalEnergy2D(Equation):
    """
        比内能变化率 书(5-26)
    """

    def __init__(self, dest, sources, alpha=1, beta=1):
        self.alpha = alpha
        self.beta = beta
        super(InternalEnergy2D, self).__init__(dest, sources)

    def initialize(self, d_ae, d_idx, d_v00, d_v01, d_v10, d_v11, d_s00, d_s01, d_s11, d_rho):
        # 应变率张量变化率 书(5-29)
        eps00 = d_v00[d_idx]
        eps01 = 0.5 * (d_v01[d_idx] + d_v10[d_idx])
        eps10 = eps01
        eps11 = d_v11[d_idx]

        # 5-26 尾项
        d_ae[d_idx] = (1 / d_rho[d_idx]) * (eps00 * d_s00[d_idx] + eps01 * d_s01[d_idx] +
                                            eps10 * d_s01[d_idx] + eps11 * d_s11[d_idx])

    def loop(self, d_idx, d_ae, DWIJ, s_m, s_idx, d_rho, s_rho, d_p, s_p, d_u, d_v, s_u, s_v):
        # 5-26 求和项，没加人工粘度项
        d_ae[d_idx] += 0.5 * (s_m[s_idx] *
                              ((d_p[d_idx] + s_p[s_idx]) / (d_rho[d_idx] * s_rho[s_idx])) *
                              ((d_u[d_idx] - s_u[s_idx]) * DWIJ[0] +
                               (d_v[d_idx] - s_v[s_idx]) * DWIJ[1]
                               )
                              )


class ContinuityEquation2D(Equation):
    """
        连续性方程
        计算密度变化率
    """

    def initialize(self, d_idx, d_arho):
        d_arho[d_idx] = 0.0

    def loop(self, d_idx, d_arho, s_idx, s_m, DWIJ, VIJ):
        vijdotdwij = DWIJ[0] * VIJ[0] + DWIJ[1] * VIJ[1]
        d_arho[d_idx] += s_m[s_idx] * vijdotdwij


class VelocityGradient2D(Equation):
    r"""
        计算速度对位置的梯度
    """

    def initialize(self, d_idx, d_v00, d_v01, d_v10, d_v11):
        d_v00[d_idx] = 0.0
        d_v01[d_idx] = 0.0
        d_v10[d_idx] = 0.0
        d_v11[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_m, s_rho, d_v00, d_v01, d_v10, d_v11, DWIJ, VIJ):
        tmp = s_m[s_idx] / s_rho[s_idx]

        d_v00[d_idx] += tmp * -VIJ[0] * DWIJ[0]
        d_v01[d_idx] += tmp * -VIJ[0] * DWIJ[1]

        d_v10[d_idx] += tmp * -VIJ[1] * DWIJ[0]
        d_v11[d_idx] += tmp * -VIJ[1] * DWIJ[1]


class RadialReturn2D(Equation):
    """ 应力退回算法 书-P155 """

    def __init__(self, dest, sources):

        super(RadialReturn2D, self).__init__(dest, sources)

    def loop(self, d_idx, d_Y, d_dam, d_s00, d_s01, d_s11, d_G):
        J2_ = (d_s00[d_idx] ** 2 + d_s11[d_idx] ** 2) / 2 + d_s01[d_idx] ** 2
        J2 = (3 * J2_ / 2) ** 0.5
        if J2 < d_Y[d_idx]:
            pass
        else:
            d_dam[d_idx] += ((J2 - d_Y[d_idx]) / (3 * d_G[d_idx]))


class HookesDeviatoricStressRate2D(Equation):
    """ 剪应力变化率 """

    def __init__(self, dest, sources, shear_mod):
        self.shear_mod = float(shear_mod)
        super(HookesDeviatoricStressRate2D, self).__init__(dest, sources)

    def initialize(self, d_idx, d_as00, d_as01, d_as11):
        d_as00[d_idx] = 0.0
        d_as01[d_idx] = 0.0

        d_as11[d_idx] = 0.0

    def loop(self, d_idx, d_G, d_D,
             d_s00, d_s01, d_s11,
             d_v00, d_v01,
             d_v10, d_v11,
             d_as00, d_as01,
             d_as11):
        if d_D[d_idx] != 1:
            v00 = d_v00[d_idx]
            v01 = d_v01[d_idx]

            v10 = d_v10[d_idx]
            v11 = d_v11[d_idx]

            s00 = d_s00[d_idx]
            s01 = d_s01[d_idx]

            s10 = d_s01[d_idx]
            s11 = d_s11[d_idx]

            # strain rate tensor is symmetric
            eps00 = v00
            eps01 = 0.5 * (v01 + v10)

            eps11 = v11

            # rotation tensor is asymmetric
            omega00 = 0.0
            omega01 = 0.5 * (v01 - v10)

            omega10 = -omega01
            omega11 = 0.0

            # tmp = 2.0*self.shear_mod
            tmp = 2.0 * d_G[d_idx]
            trace = 1.0 / 3.0 * (eps00 + eps11)

            # S_00
            d_as00[d_idx] = tmp * (eps00 - trace) + \
                            (s00 * omega00 + s01 * omega01) + \
                            (s00 * omega00 + s10 * omega01)

            # S_01
            d_as01[d_idx] = tmp * (eps01) + \
                            (s00 * omega10 + s01 * omega11) + \
                            (s01 * omega00 + s11 * omega01)

            # S_11
            d_as11[d_idx] = tmp * (eps11 - trace) + \
                            (s10 * omega10 + s11 * omega11) + \
                            (s01 * omega10 + s11 * omega11)


class MomentumEquation2D(Equation):
    """ 动量方程 不带应力张量 """

    def __init__(self, dest, sources, c0,
                 alpha=1.0, beta=1.0, gx=0.0, gy=0.0, gz=0.0,
                 tensile_correction=False):

        self.alpha = alpha
        self.beta = beta
        self.gx = gx
        self.gy = gy
        self.gz = gz
        self.c0 = c0

        self.tensile_correction = tensile_correction

        super(MomentumEquation2D, self).__init__(dest, sources)

    def initialize(self, d_idx, d_au, d_av, d_dt_cfl):
        d_au[d_idx] = 0.0
        d_av[d_idx] = 0.0
        d_dt_cfl[d_idx] = 0.0

    def loop(self, d_idx, s_idx, d_rho, d_cs,
             d_p, d_au, d_av, s_m,
             s_rho, s_cs, s_p, VIJ,
             XIJ, HIJ, R2IJ, RHOIJ1, EPS, d_D, s_D,
             DWIJ, WIJ, WDP, d_dt_cfl):

        d_s = s_D[s_idx]
        d_d = d_D[d_idx]

        rhoi21 = 1.0 / (d_rho[d_idx] * d_rho[d_idx])
        rhoj21 = 1.0 / (s_rho[s_idx] * s_rho[s_idx])

        vijdotxij = VIJ[0] * XIJ[0] + VIJ[1] * XIJ[1] + VIJ[2] * XIJ[2]

        piij = 0.0
        if vijdotxij < 0:
            cij = 0.5 * (d_cs[d_idx] + s_cs[s_idx])

            muij = (HIJ * vijdotxij) / (R2IJ + EPS)

            piij = -self.alpha * cij * muij + self.beta * muij * muij
            piij = piij * RHOIJ1

        # compute the CFL time step factor
        _dt_cfl = 0.0
        if R2IJ > 1e-12:
            _dt_cfl = abs(HIJ * vijdotxij / R2IJ) + self.c0
            d_dt_cfl[d_idx] = max(_dt_cfl, d_dt_cfl[d_idx])

        tmpi = d_p[d_idx] * rhoi21
        tmpj = s_p[s_idx] * rhoj21
        # if d_p[d_idx]<0:
        #     tmpi*=(1-d_d)
        # if s_p[s_idx]<0:
        #     tmpj*=(1-d_s)

        fij = WIJ / WDP
        Ri = 0.0
        Rj = 0.0

        # tensile instability correction
        if self.tensile_correction:
            fij = fij * fij
            fij = fij * fij

            if d_p[d_idx] > 0:
                Ri = 0.01 * tmpi
            else:
                Ri = 0.2 * abs(tmpi)

            if s_p[s_idx] > 0:
                Rj = 0.01 * tmpj
            else:
                Rj = 0.2 * abs(tmpj)

        # gradient and correction terms
        tmp = (tmpi + tmpj) + (Ri + Rj) * fij

        d_au[d_idx] += -s_m[s_idx] * (tmp + piij) * DWIJ[0]
        d_av[d_idx] += -s_m[s_idx] * (tmp + piij) * DWIJ[1]

    def post_loop(self, d_idx, d_au, d_av, d_dt_force):
        d_au[d_idx] += self.gx
        d_av[d_idx] += self.gy

        acc2 = (d_au[d_idx] * d_au[d_idx] +
                d_av[d_idx] * d_av[d_idx])

        # store the square of the max acceleration
        d_dt_force[d_idx] = acc2


class MomentumEquationWithStress2D(Equation):
    """ 应力张量的动量方程 """

    def __init__(self, dest, sources, wdeltap=-1, n=1):
        self.wdeltap = wdeltap
        self.n = n
        self.with_correction = True
        if wdeltap < 0:
            self.with_correction = False
        super(MomentumEquationWithStress2D, self).__init__(dest, sources)

    def initialize(self, d_idx, d_au, d_av):
        d_au[d_idx] = 0.0
        d_av[d_idx] = 0.0

    def loop(self, d_idx, s_idx, d_rho, s_rho, s_m,
             d_s00, d_s01, d_s11,
             s_s00, s_s01, s_s11, d_D, s_D,
             d_au, d_av, DWIJ):
        rhoa = d_rho[d_idx]
        rhob = s_rho[s_idx]

        rhoa21 = 1. / (rhoa * rhoa)
        rhob21 = 1. / (rhob * rhob)

        drate = 1 - d_D[d_idx]  ##########
        srate = 1 - s_D[s_idx]  ##########

        s00a = drate * d_s00[d_idx]
        s01a = drate * d_s01[d_idx]

        s10a = drate * d_s01[d_idx]
        s11a = drate * d_s11[d_idx]

        s00b = srate * s_s00[s_idx]
        s01b = srate * s_s01[s_idx]

        s10b = srate * s_s01[s_idx]
        s11b = srate * s_s11[s_idx]

        s00a = s00a
        s00b = s00b

        s11a = s11a
        s11b = s11b

        mb = s_m[s_idx]

        d_au[d_idx] += mb * (s00a * rhoa21 + s00b * rhob21) * DWIJ[0] + \
                       mb * (s01a * rhoa21 + s01b * rhob21) * DWIJ[1]

        d_av[d_idx] += mb * (s10a * rhoa21 + s10b * rhob21) * DWIJ[0] + \
                       mb * (s11a * rhoa21 + s11b * rhob21) * DWIJ[1]

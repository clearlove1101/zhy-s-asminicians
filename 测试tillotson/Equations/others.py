"""
    其它方程
"""
# from pysph.sph.basic_equations import ContinuityEquation
from pysph.sph.equation import Equation
from pysph.sph.solid_mech.basic import MomentumEquationWithStress
import config

# from pysph.sph.wc.basic import MomentumEquation


# from pysph.sph.solid_mech.basic import HookesDeviatoricStressRate


class InternalEnergy3D(Equation):
    """
        比内能变化率 书(5-26)
    """

    def __init__(self, dest, sources, alpha=1, beta=1):
        self.alpha = alpha
        self.beta = beta
        super(InternalEnergy3D, self).__init__(dest, sources)

    def initialize(self, d_ae, d_idx, d_v00, d_v01, d_v10, d_v11, d_v02, d_v20, d_v12, d_v21, d_v22,
                   d_s00, d_s01, d_s11, d_s02, d_s12, d_s22, d_rho):
        # 应变率张量变化率 书(5-29)

        eps00 = d_v00[d_idx]
        eps01 = 0.5 * (d_v01[d_idx] + d_v10[d_idx])
        eps10 = eps01
        eps11 = d_v11[d_idx]
        eps02 = 0.5 * (d_v02[d_idx] + d_v20[d_idx])
        eps20 = eps02
        eps12 = 0.5 * (d_v12[d_idx] + d_v21[d_idx])
        eps21 = eps12
        eps22 = d_v22[d_idx]


        # 5-26 尾项
        d_ae[d_idx] = (1 / d_rho[d_idx]) \
        * (eps00 * d_s00[d_idx] + eps01 * d_s01[d_idx] + eps02 * d_s02[d_idx] +
        eps10 * d_s01[d_idx] + eps11 * d_s11[d_idx] + eps12 * d_s12[d_idx] +
           eps20 * d_s02[d_idx] + eps21 * d_s12[d_idx] + eps22 * d_s22[d_idx])   # 多数项是相同的，可代码优化

    def loop(self, d_idx, d_ae, DWIJ, s_m, s_idx, d_rho, s_rho, d_p, s_p, d_u, d_v, d_w, s_u, s_v, s_w):
        # 5-26 求和项，没加人工粘度项
        d_ae[d_idx] += 0.5 * (s_m[s_idx] *
                              ((d_p[d_idx] + s_p[s_idx]) / (d_rho[d_idx] * s_rho[s_idx])) *
                              ((d_u[d_idx] - s_u[s_idx]) * DWIJ[0] +
                               (d_v[d_idx] - s_v[s_idx]) * DWIJ[1] +
                               (d_w[d_idx] - s_w[s_idx]) * DWIJ[2]
                               )
                              )


class ContinuityEquation3D(Equation):
    """
        连续性方程
        计算密度变化率
    """

    def initialize(self, d_idx, d_arho):
        d_arho[d_idx] = 0.0

    def loop(self, d_idx, d_arho, s_idx, s_m, DWIJ, VIJ):
        vijdotdwij = DWIJ[0] * VIJ[0] + DWIJ[1] * VIJ[1] + DWIJ[2] * VIJ[2]
        d_arho[d_idx] += s_m[s_idx] * vijdotdwij


class VelocityGradient3D(Equation):
    r"""
        计算速度对位置的梯度
    """

    def initialize(self, d_idx, d_v00, d_v01, d_v10, d_v11):
        d_v00[d_idx] = 0.0
        d_v01[d_idx] = 0.0
        d_v10[d_idx] = 0.0
        d_v11[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_m, s_rho,
             d_v00, d_v01, d_v10, d_v11, d_v02, d_v12, d_v20, d_v21, d_v22,
             DWIJ, VIJ):
        tmp = s_m[s_idx] / s_rho[s_idx]

        d_v00[d_idx] += tmp * -VIJ[0] * DWIJ[0]
        d_v01[d_idx] += tmp * -VIJ[0] * DWIJ[1]
        d_v02[d_idx] += tmp * -VIJ[0] * DWIJ[2]

        d_v10[d_idx] += tmp * -VIJ[1] * DWIJ[0]
        d_v11[d_idx] += tmp * -VIJ[1] * DWIJ[1]
        d_v12[d_idx] += tmp * -VIJ[1] * DWIJ[2]

        d_v20[d_idx] += tmp * -VIJ[2] * DWIJ[0]
        d_v21[d_idx] += tmp * -VIJ[2] * DWIJ[1]
        d_v22[d_idx] += tmp * -VIJ[2] * DWIJ[2]


class RadialReturn3D(Equation):
    """ 应力退回算法 书-P155 """

    def __init__(self, dest, sources):

        super(RadialReturn3D, self).__init__(dest, sources)

    def loop(self, d_idx, d_Y, d_dam, d_s00, d_s01, d_s11, d_s02, d_s12, d_s22, d_G):

        J2_ = (d_s00[d_idx] ** 2 + d_s11[d_idx] ** 2 + d_s22[d_idx] ** 2) / 2 \
              + d_s01[d_idx] ** 2 + d_s12[d_idx] ** 2 + d_s02[d_idx] ** 2
        J2 = (3 * J2_) ** 0.5

        if J2 < d_Y[d_idx]:
            pass
        else:
            d_dam[d_idx] += ((J2 - d_Y[d_idx]) / (3 * d_G[d_idx]))
            ratio = d_Y[d_idx] / J2

            d_s00[d_idx] *= ratio
            d_s01[d_idx] *= ratio
            d_s02[d_idx] *= ratio
            d_s11[d_idx] *= ratio
            d_s12[d_idx] *= ratio
            d_s22[d_idx] *= ratio


class HookesDeviatoricStressRate3D(Equation):
    r""" **Rate of change of stress **

    .. math::
        \frac{dS^{ij}}{dt} = 2\mu\left(\epsilon^{ij} - \frac{1}{3}\delta^{ij}
        \epsilon^{ij}\right) + S^{ik}\Omega^{jk} + \Omega^{ik}S^{kj}

    where

    .. math::

        \epsilon^{ij} = \frac{1}{2}\left(\frac{\partial v^i}{\partial x^j} +
        \frac{\partial v^j}{\partial x^i}\right)\\

        \Omega^{ij} = \frac{1}{2}\left(\frac{\partial v^i}{\partial x^j} -
           \frac{\partial v^j}{\partial x^i} \right)

    """
    def __init__(self, dest, sources, shear_mod):
        r"""
        Parameters
        ----------
        shear_mod : float
            shear modulus (:math:`\mu`)
        """

        self.shear_mod = float(shear_mod)
        super(HookesDeviatoricStressRate3D, self).__init__(dest, sources)

    def initialize(self, d_idx, d_as00, d_as01, d_as02, d_as11, d_as12, d_as22):
        d_as00[d_idx] = 0.0
        d_as01[d_idx] = 0.0
        d_as02[d_idx] = 0.0

        d_as11[d_idx] = 0.0
        d_as12[d_idx] = 0.0

        d_as22[d_idx] = 0.0

    def loop(self, d_idx, d_G,
             d_s00, d_s01, d_s02, d_s11, d_s12, d_s22,
             d_v00, d_v01, d_v02,
             d_v10, d_v11, d_v12,
             d_v20, d_v21, d_v22,
             d_as00, d_as01, d_as02,
             d_as11, d_as12,
             d_as22):

        v00 = d_v00[d_idx]
        v01 = d_v01[d_idx]
        v02 = d_v02[d_idx]

        v10 = d_v10[d_idx]
        v11 = d_v11[d_idx]
        v12 = d_v12[d_idx]

        v20 = d_v20[d_idx]
        v21 = d_v21[d_idx]
        v22 = d_v22[d_idx]

        s00 = d_s00[d_idx]
        s01 = d_s01[d_idx]
        s02 = d_s02[d_idx]

        s10 = d_s01[d_idx]
        s11 = d_s11[d_idx]
        s12 = d_s12[d_idx]

        s20 = d_s02[d_idx]
        s21 = d_s12[d_idx]
        s22 = d_s22[d_idx]

        # strain rate tensor is symmetric
        eps00 = v00
        eps01 = 0.5 * (v01 + v10)
        eps02 = 0.5 * (v02 + v20)

        eps10 = eps01
        eps11 = v11
        eps12 = 0.5 * (v12 + v21)

        eps20 = eps02
        eps21 = eps12
        eps22 = v22

        # rotation tensor is asymmetric
        omega00 = 0.0
        omega01 = 0.5 * (v01 - v10)
        omega02 = 0.5 * (v02 - v20)

        omega10 = -omega01
        omega11 = 0.0
        omega12 = 0.5 * (v12 - v21)

        omega20 = -omega02
        omega21 = -omega12
        omega22 = 0.0

        # tmp = 2.0*self.shear_mod
        tmp = 2.0 * d_G[d_idx]
        trace = 1.0/3.0 * (eps00 + eps11)

        # S_00
        d_as00[d_idx] = tmp*( eps00 - trace ) + \
                        ( s00*omega00 + s01*omega01 + s02*omega02) + \
                        ( s00*omega00 + s10*omega01 + s20*omega02)

        # S_01
        d_as01[d_idx] = tmp*(eps01) + \
                        ( s00*omega10 + s01*omega11 + s02*omega12) + \
                        ( s01*omega00 + s11*omega01 + s21*omega02)

        # S_02
        d_as02[d_idx] = tmp*eps02 + \
                        (s00*omega20 + s01*omega21 + s02*omega22) + \
                        (s02*omega00 + s12*omega01 + s22*omega02)

        # S_11
        d_as11[d_idx] = tmp*( eps11 - trace ) + \
                        (s10*omega10 + s11*omega11 + s12*omega12) + \
                        (s01*omega10 + s11*omega11 + s21*omega12)

        # S_12
        d_as12[d_idx] = tmp*eps12 + \
                        (s10*omega20 + s11*omega21 + s12*omega22) + \
                        (s02*omega10 + s12*omega11 + s22*omega12)

        # S_22
        d_as22[d_idx] = tmp*(eps22 - trace) + \
                        (s20*omega20 + s21*omega21 + s22*omega22) + \
                        (s02*omega20 + s12*omega21 + s22*omega22)


class MomentumEquation3D(Equation):
    r"""**Classic Monaghan Style Momentum Equation with Artificial Viscosity**

    .. math::

        \frac{d\mathbf{v}_{a}}{dt}=-\sum_{b}m_{b}\left(\frac{p_{b}}
        {\rho_{b}^{2}}+\frac{p_{a}}{\rho_{a}^{2}}+\Pi_{ab}\right)
        \nabla_{a}W_{ab}

    where

    .. math::

        \Pi_{ab}=\begin{cases}
        \frac{-\alpha\bar{c}_{ab}\mu_{ab}+\beta\mu_{ab}^{2}}{\bar{\rho}_{ab}} &
        \mathbf{v}_{ab}\cdot\mathbf{r}_{ab}<0;\\
        0 & \mathbf{v}_{ab}\cdot\mathbf{r}_{ab}\geq0;
        \end{cases}

    with

    .. math::

        \mu_{ab}=\frac{h\mathbf{v}_{ab}\cdot\mathbf{r}_{ab}}
        {\mathbf{r}_{ab}^{2}+\eta^{2}}\\

        \bar{c}_{ab} = \frac{c_a + c_b}{2}\\

        \bar{\rho}_{ab} = \frac{\rho_a + \rho_b}{2}

    References
    ----------
    .. [Monaghan1992] J. Monaghan, Smoothed Particle Hydrodynamics, "Annual
        Review of Astronomy and Astrophysics", 30 (1992), pp. 543-574.
    """
    def __init__(self, dest, sources, c0,
                 alpha=1.0, beta=1.0, gx=0.0, gy=0.0, gz=0.0,
                 tensile_correction=False):

        r"""
        Parameters
        ----------
        c0 : float
            reference speed of sound
        alpha : float
            produces a shear and bulk viscosity
        beta : float
            used to handle high Mach number shocks
        gx : float
            body force per unit mass along the x-axis
        gy : float
            body force per unit mass along the y-axis
        gz : float
            body force per unit mass along the z-axis
        tensilte_correction : bool
            switch for tensile instability correction (Default: False)
        """

        self.alpha = alpha
        self.beta = beta
        self.gx = gx
        self.gy = gy
        self.gz = gz
        self.c0 = c0

        self.tensile_correction = tensile_correction

        super(MomentumEquation3D, self).__init__(dest, sources)

    def initialize(self, d_idx, d_au, d_av, d_aw, d_dt_cfl):
        d_au[d_idx] = 0.0
        d_av[d_idx] = 0.0
        d_aw[d_idx] = 0.0
        d_dt_cfl[d_idx] = 0.0

    def loop(self, d_idx, s_idx, d_rho, d_cs,
             d_p, d_au, d_av, d_aw, s_m,
             s_rho, s_cs, s_p, VIJ,
             XIJ, HIJ, R2IJ, RHOIJ1, EPS,
             DWIJ, WIJ, WDP, d_dt_cfl):

        rhoi21 = 1.0/(d_rho[d_idx]*d_rho[d_idx])
        rhoj21 = 1.0/(s_rho[s_idx]*s_rho[s_idx])

        vijdotxij = VIJ[0]*XIJ[0] + VIJ[1]*XIJ[1] + VIJ[2]*XIJ[2]

        piij = 0.0
        if vijdotxij < 0:
            cij = 0.5 * (d_cs[d_idx] + s_cs[s_idx])

            muij = (HIJ * vijdotxij)/(R2IJ + EPS)

            piij = -self.alpha*cij*muij + self.beta*muij*muij
            piij = piij*RHOIJ1

        # compute the CFL time step factor
        _dt_cfl = 0.0
        if R2IJ > 1e-12:
            _dt_cfl = abs(HIJ * vijdotxij/R2IJ) + self.c0
            d_dt_cfl[d_idx] = max(_dt_cfl, d_dt_cfl[d_idx])

        tmpi = d_p[d_idx]*rhoi21
        tmpj = s_p[s_idx]*rhoj21

        fij = WIJ/WDP
        Ri = 0.0
        Rj = 0.0

        # tensile instability correction
        if self.tensile_correction:
            fij = fij*fij
            fij = fij*fij

            if d_p[d_idx] > 0:
                Ri = 0.01 * tmpi
            else:
                Ri = 0.2*abs(tmpi)

            if s_p[s_idx] > 0:
                Rj = 0.01 * tmpj
            else:
                Rj = 0.2 * abs(tmpj)

        # gradient and correction terms
        tmp = (tmpi + tmpj) + (Ri + Rj)*fij

        d_au[d_idx] += -s_m[s_idx] * (tmp + piij) * DWIJ[0]
        d_av[d_idx] += -s_m[s_idx] * (tmp + piij) * DWIJ[1]
        d_aw[d_idx] += -s_m[s_idx] * (tmp + piij) * DWIJ[2]

    def post_loop(self, d_idx, d_au, d_av, d_aw, d_dt_force):
        d_au[d_idx] += self.gx
        d_av[d_idx] += self.gy
        d_aw[d_idx] += self.gz

        acc2 = (d_au[d_idx]*d_au[d_idx] +
                d_av[d_idx]*d_av[d_idx] +
                d_aw[d_idx]*d_aw[d_idx])

        # store the square of the max acceleration
        d_dt_force[d_idx] = acc2


class MomentumEquationWithStress3D(Equation):
    r"""**Momentum Equation with Artificial Stress**

    .. math::

        \frac{D\vec{v_a}^i}{Dt} = \sum_b m_b\left(\frac{\sigma_a^{ij}}{\rho_a^2}
        +\frac{\sigma_b^{ij}}{\rho_b^2} + R_{ab}^{ij}f^n \right)\nabla_a W_{ab}

    where

    .. math::

        f_{ab} = \frac{W(r_{ab})}{W(\Delta p)}\\

        R_{ab}^{ij} = R_{a}^{ij} + R_{b}^{ij}
    """
    def __init__(self, dest, sources, wdeltap=-1, n=1):
        r"""
        Parameters
        ----------
        wdeltap : float
            evaluated value of :math:`W(\Delta p)`
        n : float
            constant
        with_correction : bool
            switch for using tensile instability correction
        """

        self.wdeltap = wdeltap
        self.n = n
        self.with_correction = True
        if wdeltap < 0:
            self.with_correction = False
        super(MomentumEquationWithStress3D, self).__init__(dest, sources)

    def initialize(self, d_idx, d_au, d_av, d_aw):
        d_au[d_idx] = 0.0
        d_av[d_idx] = 0.0
        d_aw[d_idx] = 0.0

    def loop(self, d_idx, s_idx, d_rho, s_rho, s_m,
             d_s00, d_s01, d_s02, d_s11, d_s12, d_s22,
             s_s00, s_s01, s_s02, s_s11, s_s12, s_s22,d_D, s_D,
             d_au, d_av, d_aw, DWIJ):

        rhoa = d_rho[d_idx]
        rhob = s_rho[s_idx]

        rhoa21 = 1./(rhoa * rhoa)
        rhob21 = 1./(rhob * rhob)

        drate = 1 - d_D[d_idx]   ##########
        srate = 1 - s_D[s_idx]   ##########

        s00a = drate*d_s00[d_idx]
        s01a = drate*d_s01[d_idx]
        s02a = drate*d_s02[d_idx]

        s10a = drate*d_s01[d_idx]
        s11a = drate*d_s11[d_idx]
        s12a = drate*d_s12[d_idx]

        s20a = drate*d_s02[d_idx]
        s21a = drate*d_s12[d_idx]
        s22a = drate*d_s22[d_idx]

        s00b = srate*s_s00[s_idx]
        s01b = srate*s_s01[s_idx]
        s02b = srate*s_s02[s_idx]

        s10b = srate*s_s01[s_idx]
        s11b = srate*s_s11[s_idx]
        s12b = srate*s_s12[s_idx]

        s20b = srate*s_s02[s_idx]
        s21b = srate*s_s12[s_idx]
        s22b = srate*s_s22[s_idx]

        # compute accelerations
        mb = s_m[s_idx]

        d_au[d_idx] += mb * (s00a*rhoa21 + s00b*rhob21) * DWIJ[0] + \
                       mb * (s01a*rhoa21 + s01b*rhob21) * DWIJ[1] + \
                       mb * (s02a*rhoa21 + s02b*rhob21) * DWIJ[2]

        d_av[d_idx] += mb * (s10a*rhoa21 + s10b*rhob21) * DWIJ[0] + \
                       mb * (s11a*rhoa21 + s11b*rhob21) * DWIJ[1] + \
                       mb * (s12a*rhoa21 + s12b*rhob21) * DWIJ[2]

        d_aw[d_idx] += mb * (s20a*rhoa21 + s20b*rhob21) * DWIJ[0] + \
                       mb * (s21a*rhoa21 + s21b*rhob21) * DWIJ[1] + \
                       mb * (s22a*rhoa21 + s22b*rhob21) * DWIJ[2]

from pysph.sph.integrator_step import IntegratorStep


class EulerStep(IntegratorStep):
    """Fast but inaccurate integrator. Use this for testing"""
    def stage1(self, d_idx, d_u, d_v, d_w, d_au, d_av, d_aw, d_x, d_y,
                  d_z, d_rho, d_arho, dt):
        d_u[d_idx] += dt*d_au[d_idx]
        d_v[d_idx] += dt*d_av[d_idx]

        d_x[d_idx] += dt*d_u[d_idx]
        d_y[d_idx] += dt*d_v[d_idx]

        d_rho[d_idx] += dt*d_arho[d_idx]

        """ 这里没写的物理量要补充上 """

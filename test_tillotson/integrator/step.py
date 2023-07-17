from pysph.sph.integrator_step import IntegratorStep


class EulerStep3D(IntegratorStep):

    def stage1(self, d_idx, d_u, d_v, d_w, d_au, d_av, d_aw, d_x, d_y, d_z,
               d_rho, d_arho, d_e, d_ae,
               d_s00, d_s01, d_s11, d_s02, d_s12, d_s22,
               d_as00, d_as01, d_as11, d_as02, d_as12, d_as22,
               dt):
        if -2e3 < d_x[d_idx] < 2e3 and -2e3 < d_y[d_idx] < 1e3:  # 坐标超出这个范围的就不计算了
            d_u[d_idx] += dt * d_au[d_idx]
            d_v[d_idx] += dt * d_av[d_idx]
            d_w[d_idx] += dt * d_aw[d_idx]

            d_x[d_idx] += dt * d_u[d_idx]
            d_y[d_idx] += dt * d_v[d_idx]
            d_z[d_idx] += dt * d_w[d_idx]

            d_s00[d_idx] += dt * d_as00[d_idx]
            d_s01[d_idx] += dt * d_as01[d_idx]
            d_s11[d_idx] += dt * d_as11[d_idx]
            d_s02[d_idx] += dt * d_as02[d_idx]
            d_s12[d_idx] += dt * d_as12[d_idx]
            d_s22[d_idx] += dt * d_as22[d_idx]

            d_e[d_idx] += dt * d_ae[d_idx]

            d_rho[d_idx] += dt * d_arho[d_idx]
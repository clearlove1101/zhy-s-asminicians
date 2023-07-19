import numpy as np
from pysph.base.utils import get_particle_array_wcsph
from pysph.solver.application import Application
from pysph.solver.solver import Solver
from pysph.sph.equation import Group
from pysph.base.kernels import CubicSpline
from pysph.sph.integrator import EulerIntegrator
from integrator import step
from pysph.tools import geometry as G
from Equations import EOS
import config
from visualization import visual


class App(Application):

    def create_particles(self):

        x, y = G.get_3d_sphere(dx=dx_asteroid, center=center_asteroid, r=r_asteroid)
        asteroid = get_particle_array_wcsph(
            name=name_asteroid, x=x, y=y, h=h_asteroid, m=rho_asteroid * dx_asteroid ** 2, rho=rho_asteroid,
            as00=as00_asteroid, as01=as01_asteroid, as11=as11_asteroid,
            s00=s00_asteroid, s01=s01_asteroid, s11=s11_asteroid,
            G=G_asteroid, T=T_asteroid, Y=Y_asteroid,
            dam=dam_asteroid, e=e_asteroid, ae=ae_asteroid,
            v00=v00_asteroid, v01=v01_asteroid, v10=v10_asteroid, v11=v11_asteroid,
            D=D_asteroid, cwij=cwij_asteroid
        )

        # v
        x, y = G.get_3d_sphere(dx=dx_impactor, center=center_impactor, r=r_impactor)
        impactor = get_particle_array_wcsph(
            name=name_impactor, x=x, y=y, h=h_impactor, m=rho_impactor * dx_impactor ** 2, rho=rho_impactor,
            as00=as00_impactor, as01=as01_impactor, as11=as11_impactor,
            s00=s00_impactor, s01=s01_impactor, s11=s11_impactor,
            G=G_impactor, T=T_impactor, Y=Y_impactor,
            dam=dam_impactor, e=e_impactor, ae=ae_impactor,
            v00=v00_impactor, v01=v01_impactor, v10=v10_impactor, v11=v11_impactor,
            D=D_impactor, cwij=cwij_impactor,
            v=-100
        )

        return [asteroid, impactor]

    def create_equations(self):

        equations = [Group(equations=[
            #EOS.TillotsonEOS(dest=name_impactor, sources=None, rho0=rho_impactor, c0=1.5e3, gamma=7),
            #EOS.TillotsonEOS(dest=name_asteroid, sources=None, rho0=rho_asteroid, c0=1.5e3, gamma=7),

            EOS.TaitEOS(dest=name_impactor, sources=None, rho0=rho_impactor, c0=1.5e3, gamma=7),
            EOS.TaitEOS(dest=name_asteroid, sources=None, rho0=name_asteroid, c0=1.5e3, gamma=7),
        ])]
        return equations

    def create_solver(self):

        integrator = EulerIntegrator(asteroid=step.EulerStep3D(), impactor=step.EulerStep3D())

        solver = Solver(dim=3, kernel=CubicSpline(dim=3), tf=0.5, dt=1e-7, integrator=integrator)  # 参数还得细调
        return solver


app = App()

app.run()

visual.show_all_process(r"./main_output", s=1)


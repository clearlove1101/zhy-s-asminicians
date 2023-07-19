# coding="utf-8"
from pysph.solver.application import Application
from pysph.solver.solver import Solver
from pysph.sph.equation import Group
from pysph.base.kernels import CubicSpline
from pysph.base.utils import get_particle_array_wcsph
from Equations import EOS, others, solid
from pysph.sph.integrator import EulerIntegrator
from integrator import step
from pysph.tools import geometry
from config import *
from visualization import visual
from tools import tools


class App(Application):

    # 整个仿真的类
    def create_particles(self):
        # 初始化粒子

        # 构建表示小行星的粒子
        # 球形
        x, y = geometry.get_2d_circle(dx=dx_asteroid, center=center_asteroid, r=r_asteroid)

        asteroid = get_particle_array_wcsph(name=name_asteroid, x=x, y=y,
                                            h=h_asteroid, m=rho_asteroid * dx_asteroid ** 2, rho=rho_asteroid,
                                            G=G_asteroid, T=T_asteroid, Y=Y_asteroid,
                                            as00=as00_asteroid, as01=as01_asteroid, as11=as11_asteroid,
                                            s00=s00_asteroid, s01=s01_asteroid, s11=s11_asteroid,
                                            dam=dam_asteroid, e=e_asteroid, ae=ae_asteroid,
                                            v00=v00_asteroid, v01=v01_asteroid, v10=v10_asteroid, v11=v11_asteroid,
                                            D=D_asteroid, cwij=cwij_asteroid)

        # 构建表示撞击体的粒子
        x, y = geometry.get_2d_circle(dx=dx_impactor, center=center_impactor, r=r_impactor)
        u, v = tools.get_velocity2D(V_impactor, angle_impactor)
        impactor = get_particle_array_wcsph(name=name_impactor, x=x, y=y, u=u, v=v,
                                            h=h_impactor, rho=rho_asteroid, m=rho_impactor * dx_asteroid ** 2,
                                            G=G_impactor, T=T_impactor, Y=Y_impactor,
                                            as00=as00_impactor, as01=as01_impactor, as11=as11_impactor,
                                            s00=s00_impactor, s01=s01_impactor, s11=s11_impactor,
                                            dam=dam_impactor, e=e_impactor, ae=ae_impactor,
                                            v00=v00_impactor, v01=v01_impactor, v10=v10_impactor, v11=v11_impactor,
                                            D=D_impactor, cwij=cwij_impactor)
        asteroid.add_output_arrays(["T", "s01", "v01", "as01", "av", "dam", "e", "D"])

        return [asteroid, impactor]

    def create_equations(self):
        # 给粒子加上公式，公式会约束粒子的行为

        equations = [Group(equations=[

            # EOS.TaitEOS(dest=name_impactor, sources=None, rho0=rho_impactor, c0=tillotson_c0, gamma=7),
            # EOS.TaitEOS(dest=name_asteroid, sources=None, rho0=rho_asteroid, c0=tillotson_c0, gamma=7),
            # #
            EOS.TillotsonEOS_new(dest=name_impactor, sources=None),
            EOS.TillotsonEOS_new(dest=name_asteroid, sources=None),

            others.InternalEnergy2D(dest=name_impactor, sources=[name_asteroid, name_impactor]),
            others.InternalEnergy2D(dest=name_asteroid, sources=[name_asteroid, name_impactor]),

            others.ContinuityEquation2D(dest=name_impactor, sources=[name_asteroid, name_impactor]),
            others.ContinuityEquation2D(dest=name_asteroid, sources=[name_asteroid, name_impactor]),

            others.VelocityGradient2D(dest=name_impactor, sources=[name_impactor]),
            others.VelocityGradient2D(dest=name_asteroid, sources=[name_asteroid]),

            solid.ROCK(dest=name_impactor, sources=None),
            solid.ROCK(dest=name_asteroid, sources=None),

            others.RadialReturn2D(dest=name_impactor, sources=[name_impactor]),
            others.RadialReturn2D(dest=name_asteroid, sources=[name_asteroid]),

            others.HookesDeviatoricStressRate2D(dest=name_impactor, sources=[name_impactor], shear_mod=5.6e9),
            others.HookesDeviatoricStressRate2D(dest=name_asteroid, sources=[name_asteroid], shear_mod=5.6e9),

            others.MomentumEquation2D(dest=name_impactor, sources=[name_asteroid, name_impactor], c0=tillotson_c0),
            others.MomentumEquation2D(dest=name_asteroid, sources=[name_asteroid, name_impactor], c0=tillotson_c0),

            others.MomentumEquationWithStress2D(dest=name_impactor, sources=[name_impactor]),
            others.MomentumEquationWithStress2D(dest=name_asteroid, sources=[name_asteroid]),

            #     1. InternalEnergy（比内能）
            #     2. Tillotson代替TaitEOS（别忘了计算温度）
            #     3. ROCK （适用于岩石的本构方程）
            #     4. VelocityGradient（库里有，速度梯度）
            #     5. HookesDeviatoricStressRate （剪应力变化率）
            #     6. RadialReturn （应力退回算法）
            #     7. MomentumEquationWithStress （动量方程应力张量形式）
            #     8. 孔隙度模型
            #     9. 自重力模型
        ])]
        return equations

    def create_solver(self):
        # 计算积分器（比如怎么通过当前粒子速度计算下一时刻粒子位置）

        # 欧拉法计算量最小，但精度不高（其实也还可以）
        integrator = EulerIntegrator(asteroid=step.EulerStep(), impactor=step.EulerStep())

        # CubicSpline是计算量最小的核函数（加权求和的权重）
        # tf是运行时间
        # dt是每两步间隔时间
        solver = Solver(dim=2, kernel=CubicSpline(dim=2), tf=1, dt=5e-5, integrator=integrator)  # 参数还得细调
        return solver


# 实例化类
app = App()

# 运行，运行的同时会自动保存结果
app.run()

visual.show_all_process(
    r"./main_output", s=1, xlimit=[-2 * r_asteroid, 2 * r_asteroid], ylimit=[-2 * r_asteroid, 2 * r_asteroid]
)

"""
    可视化工具，基于matplotlab
"""
import numpy as np
import matplotlib.pyplot as plt
import glob


def load_npz(path: str) -> dict:
    """
        load particles from .npz file
    :param path: .npy path
    :return: a dict for all particles
    """
    data = np.load(path, allow_pickle=True)
    particles = data['particles'].item()
    return particles


def show_one_image(
        path: str, particle_names: list = None, property_name: str = 'rho',
        gba: list = None, xlimit: list = None, ylimit: list = None, s: int = 10, immediately_show: bool = True
):
    """
        show particles of one stage from .npz file
    :param path:  .npy file path
    :param particle_names:  要可视化哪些粒子，写名字
    :param property_name:   用颜色表示粒子的哪个属性
    :param gba:             颜色表示rgba的gba
    :param xlimit:          散点图只表示x在这个范围内的粒子
    :param ylimit:          y范围
    :param s:               点大小
    :param immediately_show: 是否立刻生成可视化图像
    :return:  None
    """

    # default gba is [0.2, 0.2, 1]
    if gba is None:
        gba = [0.2, 0.2, 1]

    # load particles
    particles = load_npz(path)

    # default to show all particles
    if particle_names is None:
        particle_names = particles.keys()

    # gather all particles to show
    particles = [particles[particle]['arrays'] for particle in particle_names]

    # xyz
    x = np.concatenate([particle['x'] for particle in particles])
    y = np.concatenate([particle['y'] for particle in particles])
    z = np.concatenate([particle['z'] for particle in particles])

    # calculate color by particle property
    property = np.concatenate([particle[property_name] for particle in particles])
    # unitization
    min_p = np.min(property)
    max_p = np.max(property)
    property = (property - min_p) / (max_p - min_p + min((min_p) * 1e-10, 1e-20))
    # particle number
    num = property.shape[0]
    # rgba
    colors = np.concatenate([
        property.reshape([num, 1]),
        np.full([num, 1], gba[0]),
        np.full([num, 1], gba[1]),
        np.full([num, 1], gba[2])
    ], axis=-1)

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.scatter(x, y, z, c=colors, marker='o', s=s)
    # ax1.ylim(ylimit)

    # plt.ylim(ylimit)
    if immediately_show:
        plt.show()


def show_particles(
        objs: list, property_name: str = 'rho',
        gba: list = None, xlimit: list = None, ylimit: list = None, s: int = 10, immediately_show: bool = True
):
    """
            show particles from a list of particles
        :param objs:            粒子们
        :param property_name:   用颜色表示粒子的哪个属性
        :param gba:             颜色表示rgba的gba
        :param xlimit:          散点图只表示x在这个范围内的粒子
        :param ylimit:          y范围
        :param s:               点大小
        :param immediately_show: 是否立刻生成可视化图像
        :return:  None
        """

    # default gba is [0.2, 0.2, 1]
    if gba is None:
        gba = [0.2, 0.2, 1]

    # gather all particles to show
    particles = objs

    # xy
    x = np.concatenate([particle.x for particle in particles])
    y = np.concatenate([particle.y for particle in particles])
    z = np.concatenate([particle.z for particle in particles])

    # calculate color by particle property
    """ 需要其它属性的话得手动加 """
    if property_name == 'rho':
        property = np.concatenate([particle.rho for particle in particles])
    elif property_name == 'p':
        property = np.concatenate([particle.p for particle in particles])
    elif property_name == 'u':
        property = np.concatenate([particle.u for particle in particles])

    # unitization
    min_p = np.min(property)
    max_p = np.max(property)
    property = (property - min_p) / (max_p - min_p + min((min_p) * 1e-10, 1e-20))
    # particle number
    num = property.shape[0]
    # rgba
    colors = np.concatenate([
        property.reshape([num, 1]),
        np.full([num, 1], gba[0]),
        np.full([num, 1], gba[1]),
        np.full([num, 1], gba[2])
    ], axis=-1)

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.scatter(x, y, z, c=colors, marker='o', s=s)
    if immediately_show:
        plt.show()


def demo():
    from pysph.base.utils import get_particle_array_wcsph
    from pysph.tools import geometry as G

    # create particles
    x, y = G.get_2d_block(dx=1, length=100, height=100, center=[10, 10])
    grd = get_particle_array_wcsph(name='grd', x=x, y=y, h=1.2, m=1000, rho=1000)

    # visualize the particles
    show_particles([grd], s=1)


def show_all_process(
        dir_path: str, num: list = None, particle_names: list = None, property_name: str = 'rho',
        gba: list = None, xlimit: list = None, ylimit: list = None, s: int = 10
):
    """
        从文件夹里输出一些时刻的图像
    :param dir_path:  输出文件夹
    :param num:       图像输出几排几个
    :param particle_names:  下面的同def show_one_image()
    :param property_name:
    :param gba:
    :param xlimit:
    :param ylimit:
    :param s:
    :return: None
    """
    if num is None:
        num = [2, 3]

    lst = glob.glob(dir_path + "/*.npz")
    lst.sort(key=lambda x: int(x.split('.')[-2].split('_')[-1]))
    rate = len(lst) // (num[0] * num[1])
    rate = 1 if rate == 0 else rate

    sub_num = len(lst) // rate
    to_show = []
    for i in range(sub_num - 1):
        to_show.append(lst[i * rate])
    to_show.append(lst[-1])

    for idx, pth in enumerate(to_show):
        show_one_image(pth, particle_names, property_name, gba, xlimit, ylimit, s, immediately_show=True)

    print("plots from {}".format(dir_path))
    # plt.show()

# demo()
# show_all_process(r"../../asteroid/Dec4_22_output")
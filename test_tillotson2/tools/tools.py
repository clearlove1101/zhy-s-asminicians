
"""
    给粒子加自转角速度等函数
"""
import numpy as np


def add_omg(x: np.array, y: np.array, z: np.array, omg: list, center: list = None) -> [np.array, np.array, np.array]:
    """
        计算每个粒子由于全局的角速度增加的速度
    :param x: 粒子横坐标
    :param y: 粒子y坐标
    :param w: 全局角速度
    :param center: 瞬时旋转中心
    :return: 粒子增加的x、y方向速度
    """
    if center is None:
        center = [np.mean(x), np.mean(y), np.mean(z)]

    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]

    u = 0
    v = dz * omg[0]
    w = -dy * omg[0]

    u += -dz * omg[1]
    w += dx * omg[1]

    u += dy * omg[2]
    v += -dx * omg[2]

    return u, v, w



def gravity(x, y, z, m, x0, y0, z0):
    """
        计算某点引力场
    :param x:  所有粒子横坐标
    :param y:  所有粒子纵坐标
    :param m:  所有粒子质量
    :param x0: 所求引力场点横坐标
    :param y0: 所求引力场点纵坐标
    :return: 引力场两个分量
    """
    G = 6.67e-11
    dx = x - x0
    dy = y - y0
    dz = z - z0
    r = (dx ** 2 + dy ** 2 + dz) ** 0.5

    F = G * m / r ** 3
    return F * dx, F * dy, F * dz


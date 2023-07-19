"""
    给粒子加自转角速度等函数
"""
import numpy as np
import math


def add_omg(x: np.array, y: np.array, w: int, center: list = None) -> [np.array, np.array]:
    """
        计算每个粒子由于全局的角速度增加的速度
    :param x: 粒子横坐标
    :param y: 粒子y坐标
    :param w: 全局角速度
    :param center: 瞬时旋转中心
    :return: 粒子增加的x、y方向速度
    """
    if center is None:
        center = [np.mean(x), np.mean(y)]

    dx = x - center[0]
    dy = y - center[1]

    return -w * dy, w * dx


def gravity(x, y, m, x0, y0):
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
    r = (dx ** 2 + dy ** 2) ** 0.5

    F = G * m / r ** 3
    return F * dx, F * dy


def get_velocity2D(V, angle):
    u = V * math.cos(angle / 180 * math.pi)
    v = - V * math.sin(angle / 180 * math.pi)
    return u, v

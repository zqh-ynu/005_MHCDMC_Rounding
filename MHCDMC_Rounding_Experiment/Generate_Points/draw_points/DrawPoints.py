import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
import numpy as np


def draw_all_points(U_loc, A_loc, r, fname, isr=0):
    # 参数：
    # U_loc     list    [[x1, y1], ..., [xn, yn]]
    # A_loc     list    [[x1, y1], ..., [xm, ym]]
    # r         float
    # fname     str     图片存储位置
    # isr       int     是否写文件，默认为0，即否
    fig = plt.figure(figsize=(4, 4), dpi=300)
    ax = fig.add_subplot(111)
    ax.set_aspect(1)    # 设置坐标轴同比例
    # 设置坐标轴范围
    plt.xlim(-r, 2500)
    plt.ylim(-r, 2500)

    i = 0
    for a in A_loc:
        # print("a({0}, {1})".format(a[0], a[1]))
        ax.scatter(a[0], a[1], s=40, color="green", marker="1", linewidths=1, alpha=0.5)
        # radius:半径 fill:False/True,是否有填充
        # linestyle:{'-', '--', '-.', ':', '', (offset, on-off-seq), ...},线条风格
        # c = Circle((a[0], a[1]), radius=r, fill=False, linestyle='-.')
        # ax.add_patch(c)
        i += 1

    j = 0
    for u in U_loc:
        # print("u(%d, %d)" % (u[0], u[1]))
        x = u[0]
        y = u[1]
        ax.scatter(u[0], u[1], s=10, marker='*', color="blue")
        c = Circle((u[0], u[1]), radius=r, fill=False, linestyle=':')
        ax.add_patch(c)
        plt.text(x, y, str(j), fontsize=6)
        j += 1
    if isr == 1:
        plt.tight_layout()  # 去除pdf周围白边
        plt.savefig(fname)
    plt.show()

    pass

def draw_all_elements(U, A, r, fname, isr=0):
    # 参数：
    # U     list        [[x1, y1, BR1, aid], ..., [xn, yn, BRn, aid]
    # A     list        [[x1, y1, BW1, isSelect], ..., [xm, ym, BWm, isSelect]
    # r     double      半径
    # fname str         存储pdf文件路径
    # isr   int         是否写文件，默认为0，即否
    fig = plt.figure(figsize=(4, 4), dpi=300)
    ax = fig.add_subplot(111)
    ax.set_aspect(1)  # 设置坐标轴同比例
    # 设置坐标轴范围
    plt.xlim(-r, 2500)
    plt.ylim(-r, 2500)

    i = 0
    for a in A:
        x = a[0]
        y = a[1]
        isSlt = a[3]    # 是否被选择
        color="green"
        if isSlt:
            color="red"
            ax.scatter(x, y, s=40, color=color, marker="^", linewidths=1)
            plt.text(x, y, str(i), fontsize=8)
            # c = Circle((x, y), radius=r, fill=False, linestyle='-.')
            # ax.add_patch(c)
        # else:
        #     ax.scatter(x, y, s=40, color=color, marker="1", linewidths=1, alpha=0.5)
        i += 1
    j = 0
    for u in U:
        x = u[0]
        y = u[1]
        ax.scatter(x, y, s=10, marker='*', color="blue")
        # plt.text(x, y, str(j), fontsize=6)
        j += 1
    if isr == 1:
        if not os.path.exists(os.path.dirname(fname)):
            # if the demo_folder directory is not present
            # then create it.
            os.makedirs(os.path.dirname(fname))
        plt.tight_layout()  # 去除pdf周围白边

        plt.savefig(fname)
    plt.show()
    pass

# test
# from generate_points.u_random import ge_point_random
# from generate_points.a_neibor import ge_a_neibor
#
# r = 5
# U_loc = ge_point_random(10, 20)
# # print("U1\n", str(U_loc))
# A_loc = ge_a_neibor(U_loc, r)
# # print("U2\n", str(U_loc))
# draw_all_points(U_loc, A_loc, r)
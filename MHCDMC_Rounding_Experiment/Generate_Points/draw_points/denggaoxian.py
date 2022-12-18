import numpy as np
import math
import matplotlib.pyplot as plt

f = 2e9
c = 299792458
N0 = -174 + 10 * math.log10(5e5)
p = 40

print(N0)
def draw_contour(A_loc, U_loc, r, Length, fname):
    """
    将所有用户与被选择的服务器输出成图片，并画出全局的信噪比等高线图

    :param A_loc: [[x1, y1, BW1, isSelect], ..., [xm, ym, BWm, isSelect]
    :param U_loc: [[x1, y1, BR1, aid], ..., [xn, yn, BRn, aid]
    :param r: 服务器默认覆盖半径
    :param Length: 用户服务器分布区域边长
    :param fname: 图片文件保存绝对路径
    :return:
    """

    h = 100
    muban = np.linspace(-r, Length + r, int((Length + 2 * r) / 20))
    shapes = len(muban)
    # X, Y 是将分布区域的点栅格化，对应每个点的x，y坐标
    X, Y = np.meshgrid(muban, muban)

    # 计算每个点距离最近的服务器
    # 将X，Y转化为一维数组
    X = X.flatten()
    Y = Y.flatten()
    # 每个栅格点离得最近的服务器编号
    CLO_Aid = np.zeros(X.shape, dtype=np.int32)
    # 每个栅格点距离每个服务器的距离
    D = np.zeros((len(A_loc), shapes * shapes))
    for i in range(len(A_loc)):
        a = A_loc[i]
        x = a[0]
        y = a[1]
        D[i] = np.sqrt(np.power(X - x, 2) + np.power(Y - y, 2) + np.power(h, 2))
    D = D.T
    # print(D.shape)
    # print(D)
    for i in range(X.shape[0]):
        # print(type(np.argmin(D[i])))
        CLO_Aid[i] = np.argmin(D[i])  # 获取每个栅格点离得最近的服务器编号

    CLO_Aid = CLO_Aid.reshape(shapes, shapes)
    X = X.reshape(shapes, shapes)
    Y = Y.reshape(shapes, shapes)

    # 所有点到对应服务器的距离
    ALL_D = np.zeros((len(A_loc), shapes, shapes))
    # 所有点到对应服务器的损失
    ALL_L = np.zeros((len(A_loc), shapes, shapes))
    for i in range(len(A_loc)):
        a = A_loc[i]
        xa = a[0]
        ya = a[1]

        ALL_D[i] = np.sqrt(np.power(X - xa, 2) + np.power(Y - ya, 2) + np.power(h, 2))
        ALL_L[i] = 20 * np.log10(ALL_D[i] * 4 * math.pi * f / c) + 1
    # print(ALL_L)
    SINR = np.zeros((shapes, shapes))
    for i in range(shapes):
        for j in range(shapes):
            aid = CLO_Aid[i][j]
            L = ALL_L[aid][i][j]

            SINR[i][j] = p - L
            # print(SINR[i][j])
            N = math.pow(10, N0 / 10)
            # print("N1: " + str(N))
            for aa in range(len(A_loc)):
                if aa != aid:
                    L0 = ALL_L[aa][i][j]
                    N += math.pow(10, (p - L0 - 5) / 10)
            N = 10 * math.log10(N)
            # print("N2: " + str(N))
            SINR[i][j] -= N
    # print(SINR)
    levels = np.linspace(np.min(SINR), np.max(SINR), 15)

    fig, ax = plt.subplots()
    ax.set_aspect(1)  # 设置坐标轴同比例
    # 设置坐标轴范围
    plt.xlim(-r, Length + r)
    plt.ylim(-r, Length + r)

    ctf = ax.contourf(X, Y, SINR, levels=levels, alpha=0.75)
    ct = ax.contour(X, Y, SINR, levels=levels, linewidths=0.4)
    ax.clabel(ct, inline=True, fontsize=4)
    plt.colorbar(ctf)  # 添加cbar

    i = 0
    for a in A_loc:
        x = a[0]
        y = a[1]

        color = "blue"
        ax.scatter(x, y, s=40, color=color, marker="^", linewidths=1)
        # plt.text(x, y, str(i), fontsize=8)
        # else:
        #     ax.scatter(x, y, s=40, color=color, marker="1", linewidths=1, alpha=0.5)
        i += 1

    j = 0
    for u in U_loc:
        # print("u(%d, %d)" % (u[0], u[1]))
        x = u[0]
        y = u[1]
        ax.scatter(x, y, s=10, marker='*', color="red")
        plt.text(x, y, str(j), fontsize=6)
        j += 1

    plt.tight_layout()  # 去除pdf周围白边
    print(fname)
    plt.savefig(fname)

    plt.show()

    pass

# fname = "SINR.pdf"
# b = 100
# A_loc = np.array([
#     [13, 3],
#     [4, 15],
#     [10, 10],
#     [5, 15],
#     [15, 14]])
# A_loc = A_loc * b
# U_loc = []
# r = 200
# length = 20 * b
# print(A_loc)
# draw_contour(A_loc, U_loc, r, length, fname)

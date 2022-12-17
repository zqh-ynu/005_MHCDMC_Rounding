import numpy as np


def generateUserLocationRandomly(n, Length):
    """
    以随机方式生成用户坐标

    :param n: 用户数量
    :param Length: 用户分布区域边长
    :return: np.random.randint(0, Length, (n, 2))
    """
    # 按照随机方式在边长为len的区域中生成n个用户的坐标
    return np.random.randint(0, Length, (n, 2))


def generageUserLocationJ():
    """

    :return:
    """
    pass


def generateServerLocationByNeibor(U_loc, r):
    """
    根据U_loc中的用户坐标，两个用户为圆心，r为半径的圆的交点，生成服务器的坐标

    :param U_loc: 用户坐标列表[[x1, y1], [x2, y2], ..., [xn, yn]]
    :param r: 服务用户的半径
    :return: A_loc = [[x1, y1], [x2, y2], ..., [xm, ym]]
    """
    A_loc = []
    for i in range(len(U_loc)):
        # 用户1的坐标
        u1 = U_loc[i]
        A_loc.append(u1)
        # print("i", i, ":", u1)
        for j in range(i + 1, len(U_loc)):
            # 用户2的坐标
            u2 = U_loc[j]
            # print("j", j, ":", u2)
            # 用户1和用户2的距离
            d = pow(pow(u1[0] - u2[0], 2) + pow(u1[1] - u2[1], 2), 1 / 2)
            # print("d=", d)
            if d == 0:
                continue
            elif d <= 2 * r:
                # 用户1和用户2的连线中点
                u0 = [(u1[0] + u2[0]) / 2, (u1[1] + u2[1]) / 2]
                # print("uo"+str(u0))
                # 圆交点到用户1用户2连线的距离
                h = pow(pow(r, 2) - pow(d / 2, 2), 1 / 2)
                # print("h=", h)
                # 两个交点，即服务器坐标
                a1 = [0, 0]
                a2 = [0, 0]
                # 两个服务器的x坐标
                a1[0] = u0[0] - (h / d) * (u2[1] - u1[1])
                a2[0] = u0[0] + (h / d) * (u2[1] - u1[1])
                # 两个服务器的y坐标
                a1[1] = u0[1] + (h / d) * (u2[0] - u1[0])
                a2[1] = u0[1] - (h / d) * (u2[0] - u1[0])

                # print("u{0}({1},{2})与u{3}({4},{5})之间的距离为{6}，"
                #       "他们的交点坐标a1({7},{8})与a2({9},{10})"
                #       .format(i, u1[0], u1[1], j, u2[0], u2[1], d, a1[0], a1[1], a2[0], a2[1]))
                # print("u", i, "u", j, ")1:", a1)
                # print("u", i, "u", j, ")2:", a2)
                A_loc.append(a1)
                A_loc.append(a2)

    return A_loc


def generateServerLocationAverage(Length, l):
    """
    在边长为L的区域中，按照均匀分布的方式生成无人机候选位置，无人机之间相距l.从左下角坐标为(l, l)点开始部署

    :param Length: 无人机分布区域边长
    :param l: 无人机分布间隔
    :return: A_loc = [[x1, y1], [x2, y2], ..., [xm, ym]]
    """
    m = int((Length - 2 * l) / l)
    A_loc = []

    for i in range(m):
        y = (i + 1) * l
        for j in range(m):
            x = (j + 1) * l
            A_loc.append([x, y])
        x = (m + 1) * l
        A_loc.append([x, y])
    x = (m + 1) * l
    y = (m + 1) * l
    A_loc.append([x, y])
    return A_loc

# U_loc = generateUserLocationRandomly(10, 20)
#
# A_loc = generateServerLocationAverage(2000, 150)
#
# print(A_loc)

import numpy as np
import random
from generate_points.GenerateLocation import generateUserLocationRandomly


def generateUserRequest(U_loc, low_request=1, high_request=5):
    """
    为用户随机生成需求，位于区间low_request~high_request

    :param U_loc: 包含所有用户坐标的列表[[x1, y1], ..., [xn, yn]]
    :param low_request: 用户需求下界
    :param high_request: 用户需求上解
    :return: U = [[x1, y1, BR1], ..., [xn, ynBRn]]
    """

    n = len(U_loc)
    requests = np.random.randint(low_request, high_request, (n, 1))
    # 矩阵拼接，将U_loc与Requests横向拼接
    U = np.append(U_loc, requests, axis=1)
    return U


# # test:
# U_loc = generateUserLocationRandomly(10, 20)1
# print(generateUserRequest(U_loc))

def generateServerCapacity(A_loc, capacity):
    """
    为所有服务器设置相同想容量，为capacity

    :param A_loc: 包含所有服务器坐标的列表，[[x1, y1], ..., [xm, ym]]
    :param capacity: 服务器容量值
    :return: [[x1, y1, BW], ..., [xm, ym, BW]]
    """

    m = len(A_loc)

    A_loc = np.array(A_loc)
    # 初始化capacities，mx1, 初值为capacity
    capacities = np.full((m, 1), capacity)
    # 矩阵拼接，将A_loc与Capacity横向拼接
    # print(A_loc.shape)
    # print(capacities.shape)
    A = np.append(A_loc, capacities, axis=1)
    return A

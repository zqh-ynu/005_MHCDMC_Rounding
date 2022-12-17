import numpy as np
import random
from generate_points.GenerateLocation import generateUserLocationRandomly


def generateUserRequest(U_loc, low_request=1, high_request=5):
    # 参数：
    # U_loc     list        [[x1, y1], ..., [xn, yn]]
    # low_request=0
    # high_request=20

    n = len(U_loc)
    requests = np.random.randint(low_request, high_request, (n, 1))
    # 矩阵拼接，将U_loc与Requests横向拼接
    U = np.append(U_loc, requests, axis=1)
    return U


# # test:
# U_loc = generateUserLocationRandomly(10, 20)1
# print(generateUserRequest(U_loc))

def generateServerCapacity(A_loc, capacity):
    # 参数：
    # A_loc     list        [[x1, y1], ..., [xm, ym]]
    # low_capacity=20
    # high_capacity=100

    m = len(A_loc)
    # 初始化capacities，mx1, 初值为capacity
    capacities = np.full((m, 1), capacity)
    # 矩阵拼接，将A_loc与Capacity横向拼接
    A = np.append(A_loc, capacities, axis=1)
    return A

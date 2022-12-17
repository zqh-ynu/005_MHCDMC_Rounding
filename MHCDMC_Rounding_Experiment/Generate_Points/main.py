import numpy as np
import math
import matplotlib.pyplot as plt
from oneResult import readResult
import os

f = 2e9
c = 299792458
N0 = -174 + 10 * math.log10(5e5)
p = 40
h = 100

def cal_sinr(fname, n):
    # 获取项目路径
    ffname = os.path.abspath(__file__)  # 当前文件绝对路径
    fpath = os.path.dirname(ffname)  # 当前文件的父文件夹绝对路径
    fname_in = fpath + "\\data\\example\\oneInstance" + str(n) + ".txt"
    fname_out = "D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\" + fname
    A, U = readResult(fname_in, fname_out)
    A_loc = []
    AID = []
    for i in range(len(A)):
        if A[i][3] == 1:
            AID.append(i)
    for a in A:
        if a[3] == 1:
            A_loc.append(a)

    XA = np.zeros((len(A), len(U)))
    YA = np.zeros((len(A), len(U)))
    XU = np.zeros((len(A), len(U)))
    YU = np.zeros((len(A), len(U)))

    for i in range(len(A)):
        a0 = A[i]
        XA[i] += a0[0]
        YA[i] += a0[1]
    for j in range(len(U)):
        u = U[j]
        XU[0][j] += u[0]
        YU[0][j] += u[1]
    for i in range(1, len(A)):
        XU[i] += XU[0]
        YU[i] += YU[0]

    D = np.sqrt(np.power(XA - XU, 2) + np.power(YA - YU, 2) + math.pow(h, 2))
    L = 20 * np.log10(D * 4 * math.pi * f / c) + 1

    N = np.zeros((len(A), len(U)))
    for i in range(len(A)):
        for j in range(len(U)):
            n = math.pow(10, N0 / 10)
            for a in AID:
                if a != i:
                    n += math.pow(10, (p - L[a][j]) / 10)
            n = 10 * math.log10(n)
            N[i][j] = n

    SINR = p - L
    SINR -= N
    print(A_loc)
    print("-----D-----")
    print(D)
    print("-----L-----")
    print(L)
    print("-----N-----")
    print(N)
    print("-----SINR-----")
    print(SINR)
    print(SINR[2, 15])




cal_sinr("test\\MultiItTest\\n20BW50\\r1.txt", 20)

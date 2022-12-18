import os
from draw_points.DrawPoints import draw_all_elements
from draw_points.DrawPoints import draw_all_points


def readResult(infname, outfname):
    """
    从infname文件读取用户与服务器的原始数据，从outfname文件读取服务器选中数据及服务器覆盖用户的数据。
    infname文件格式：
        第1行：        m n
        第2~m+1行：    x y BW
        第m+2~m+n+1行：x y BR
    outfname文件格式：
        第1行：        cm n
        第2k+2行：     aid unum
        第2k+3行：     uid1, ..., uid_unum


    :param infname: 用户与服务器的原始数据文件
    :param outfname: 算法执行结果
    :return: A = [[x, y, BW, isSelected], ..., ], U = [[x, y, BR, servedBy], ..., ]
    """
    U=[]        # [[x, y, BR, servedBy], ..., ]
    A=[]        # [[x, y, BW, isSelected], ..., ]

    # 读取infname
    with open(infname, encoding='utf-8') as file_in:
        lines_in = file_in.readlines()
    file_in.close()
    line0 = lines_in[0][0:-1].split()
    m = int(line0[0])
    n = int(line0[1])
    for i in range(m):
        # line示例
        # 'x y BW'  x坐标，y坐标，容量BW
        line = lines_in[i + 1][0:-1].split()
        a = []
        for item in line:
            a.append(float(item))
        # 再在最后加上是否被选中的标记0为‘未选中’，1为‘选中’
        a.append(0)
        A.append(a)
    for j in range(n):
        # line示例
        # 'x y BR'  x坐标，y坐标，容量BW
        line = lines_in[j + m + 1][0:-1].split()
        u = []
        for item in line:
            u.append(float(item))
        # 再在最后加上被哪个服务器覆盖
        u.append(-1)
        U.append(u)

    # 读取outfname
    with open(outfname, encoding='utf-8') as file_out:
        lines_out = file_out.readlines()
    file_out.close()
    line0 = lines_out[0][0:-1].split()
    cm = int(line0[0])  # 被选中的服务器数量
    for k in range(cm):
        # line1 示例
        # aid unum      aid: 被选中服务器id   unum: 被该服务器服务的用户数
        line1 = lines_out[2 * k + 1][0:-1].split()
        aid = int(line1[0])
        print("a" + str(aid) + str(A[aid]))
        A[aid][3] = 1
        line2 = lines_out[2 * k + 2][0:-1].split()
        for item in line2:
            uid = int(item)
            U[uid][3] = aid
    return A, U





import os
from draw_points.DrawPoints import draw_all_elements
from draw_points.DrawPoints import draw_all_points


def oneResult(infname, outfname):
    # infname       源数据文件绝对路径
    # outfname      结果数据文件绝对路径
    U=[]
    A=[]

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
        # 'x y BW'  x坐标，y坐标，容量BW
        line = lines_in[j + m + 1][0:-1].split()
        u = []
        for item in line:
            u.append(float(item))
        # 再在最后加上被哪个服务器覆盖
        u.append(-1)
        U.append(u)

    # 读取infname
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





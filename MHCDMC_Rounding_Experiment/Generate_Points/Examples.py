import os.path

from generate_points.GenerateLocation import generateUserLocationRandomly
from generate_points.GenerateLocation import generateServerLocationByNeibor
from generate_points.GenerateLocation import generateServerLocationAverage
from generate_points.GenerateRequestAndCapacity import generateUserRequest
from generate_points.GenerateRequestAndCapacity import generateServerCapacity
from draw_points.DrawPoints import draw_all_points, draw_all_elements
from draw_points.denggaoxian import draw_contour
from PrintToFile.PrintDataToFile import print_data_to_file
from oneResult import readResult


def one_instance(n, Lenth, r, fname, Utype=1, Atype=1, l=200):
    """
    根据参数，生成问题实例

    :param n: 用户数
    :param Lenth: 用户分布区域长度
    :param r: 无人机默认覆盖半径
    :param fname: 该实例的文件绝对路径
    :param Utype: 生成用户分布的类型, 1为随机分布，2为聚集分布。默认为1
    :param Atype: 生成无人机分布的类型，1为均匀分布，2为交点分布。默认为1
    :param l: 无人机以均匀分布时的间隔距离
    :return:
    """

    if Utype == 1:
        U_loc = generateUserLocationRandomly(n, Lenth)
    elif Utype == 2:
        U_loc = generateUserLocationRandomly(n, Lenth)
    else:
        U_loc = generateUserLocationRandomly(n, Lenth)

    if Atype == 1:
        A_loc = generateServerLocationAverage(Lenth, l)
    elif Atype == 2:
        A_loc = generateServerLocationByNeibor(U_loc, r)
    else:
        A_loc = []

    U = generateUserRequest(U_loc)
    A = generateServerCapacity(A_loc, 50)

    print_data_to_file(U, A, fname)


def oneResult(fname_in, fname_out, fpath, rfile, n):
    """
    根据源数据与结果画图

    :param fname_in: 源数据
    :param fname_out: 结果数据
    :param fpath: python项目目录：D:\Myschool\graduate_school\02Graduate\Research\My paper\2_Papers\005_MHCDMC_Rounding\MHCDMC_Rounding_Experiment\Generate_Points\
    :param rfile: 结果数据目录
    :param n: 用户数
    :return: null
    """
    # 获取项目路径

    A, U = readResult(fname_in, fname_out)
    print("u5" + str(U[5]))
    print("u10" + str(U[10]))
    print("u12" + str(U[12]))
    print('m=%d, n=%d' % (len(A), len(U)))
    for a in A:
        print(a)
    print('------------------------------------')
    for u in U:
        print(u)
    print('------------------------------------')

    draw_all_elements(U, A, 500, fpath + "\\data\\example\\fig\\" + rfile[:-4] + ".pdf", 1)
    draw_all_points(U, A, 500, fpath + "\\data\\example\\fig\\oneInstance.pdf")
    A_loc = []
    for a in A:
        if a[3] == 1:
            A_loc.append(a)
    draw_contour(A_loc, U, 500, 2000, fpath + "\\data\\example\\fig\\" + rfile[:-4] + "d.pdf")


def test():
    # 获取项目路径
    ffname = os.path.abspath(__file__)  # 当前文件绝对路径
    fpath = os.path.dirname(ffname)  # 当前文件的父文件夹绝对路径
    fname = fpath + "\\data\\example\\UAV_AVEn200l200.txt"
    print(fname)
    # one_instance(200, 2000, 400, fname)

    rfile = "test\\MultiItTest\\n200BW50\\UAV_AVEn200l200_IP.txt"
    fname_in = fname
    fname_out = "D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\" + rfile
    print(fname_in)
    print(fname_out)
    oneResult(fname_in, fname_out, fpath, rfile, 200)
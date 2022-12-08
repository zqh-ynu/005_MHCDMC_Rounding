import os.path

from generate_points.GenerateLocation import generateUserLocationRandomly
from generate_points.GenerateLocation import generateServerLocationByNeibor
from generate_points.GenerateRequestAndCapacity import generateUserRequest
from generate_points.GenerateRequestAndCapacity import generateServerCapacity
from draw_points.DrawPoints import draw_all_points, draw_all_elements
from PrintToFile.PrintDataToFile import print_data_to_file
from oneResult import  oneResult

def one_instance20():
    n = 20
    len = 2000
    r = 500
    U_loc = generateUserLocationRandomly(n, len)
    A_loc = generateServerLocationByNeibor(U_loc, r)


    U = generateUserRequest(U_loc, 1, 10)
    A = generateServerCapacity(A_loc, 20)

    # 获取项目路径
    ffname = os.path.abspath(__file__)   # 当前文件绝对路径
    fpath = os.path.dirname(ffname)      # 当前文件的父文件夹绝对路径
    ppath = os.path.dirname(fpath)      # 当前文件的祖文件夹绝对路径，即项目路径

    fname = ppath + "\\data\\example\\oneInstance.txt"
    print(fname)
    print_data_to_file(U, A, fname)
    pfname = ppath + "\\data\\example\\fig\\oneInstance.pdf"
    draw_all_points(U_loc, A_loc, r, pfname)
    pass

def oneResult20(fname):
    # 获取项目路径
    ffname = os.path.abspath(__file__)  # 当前文件绝对路径
    fpath = os.path.dirname(ffname)  # 当前文件的父文件夹绝对路径

    fname_in = fpath + "\\data\\example\\oneInstance.txt"
    fname_out = "D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\" + fname
    A, U = oneResult(fname_in, fname_out)
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
    draw_all_elements(U, A, 500, fpath + "\\data\\example\\fig\\" + fname[:-4] + ".pdf", 1)
    # draw_all_points(U, A, 500, fpath + "\\data\\example\\fig\\oneInstance.pdf")

def one_instance200():
    n = 200
    len = 2000
    r = 500
    U_loc = generateUserLocationRandomly(n, len)
    A_loc = generateServerLocationByNeibor(U_loc, r)


    U = generateUserRequest(U_loc, 1, 10)
    A = generateServerCapacity(A_loc, 50)

    # 获取项目路径
    ffname = os.path.abspath(__file__)   # 当前文件绝对路径
    fpath = os.path.dirname(ffname)      # 当前文件的父文件夹绝对路径

    fname = fpath + "\\data\\example\\oneInstance200.txt"
    print(fname)
    print_data_to_file(U, A, fname)
    pfname = fpath + "\\data\\example\\fig\\oneInstance200.pdf"
    # draw_all_points(U_loc, A_loc, r, pfname)
    pass

def oneResult200():
    # 获取项目路径
    ffname = os.path.abspath(__file__)  # 当前文件绝对路径
    fpath = os.path.dirname(ffname)  # 当前文件的父文件夹绝对路径

    fname_in = fpath + "\\data\\example\\oneInstance200.txt"
    fname_out = "D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\oneresult200.txt"
    A, U = oneResult(fname_in, fname_out)
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
    draw_all_elements(U, A, 500, fpath + "\\data\\example\\fig\\oneResult200.pdf")
    # draw_all_points(U, A, 500, fpath + "\\data\\example\\fig\\oneInstance200.pdf")

# oneResult20("example1\\BW30LP.txt")
# oneResult20("example1\\BW30DSIS.txt")
# oneResult20("example1\\BW30COS.txt")

def one_instance50():
    n = 50
    len = 2000
    r = 500
    U_loc = generateUserLocationRandomly(n, len)
    A_loc = generateServerLocationByNeibor(U_loc, r)


    U = generateUserRequest(U_loc, 1, 10)
    A = generateServerCapacity(A_loc, 50)

    # 获取项目路径
    ffname = os.path.abspath(__file__)   # 当前文件绝对路径
    fpath = os.path.dirname(ffname)      # 当前文件的父文件夹绝对路径
    ppath = os.path.dirname(fpath)      # 当前文件的祖文件夹绝对路径，即项目路径

    fname = ppath + "\\data\\example\\oneInstance" + str(n) + ".txt"
    print(fname)
    print_data_to_file(U, A, fname)
    pfname = ppath + "\\data\\example\\fig\\oneInstance" + str(n) + ".pdf"
    draw_all_points(U_loc, A_loc, r, pfname)
    pass

def oneResult50(fname):
    # 获取项目路径
    ffname = os.path.abspath(__file__)  # 当前文件绝对路径
    fpath = os.path.dirname(ffname)  # 当前文件的父文件夹绝对路径

    n = 50

    fname_in = fpath + "\\data\\example\\oneInstance" + str(n) + ".txt"
    fname_out = "D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\" + fname
    A, U = oneResult(fname_in, fname_out)
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
    draw_all_elements(U, A, 500, fpath + "\\data\\example\\fig\\" + fname[:-4] + ".pdf", 1)
    draw_all_points(U, A, 500, fpath + "\\data\\example\\fig\\oneInstance.pdf")

# oneResult50("n50BW50\\rLP.txt")
# oneResult50("n50BW50\\rDSIS.txt")
# oneResult50("n50BW50\\rCOS.txt")
oneResult50("n50BW50\\rSFS.txt")
# oneResult50("n50BW50\\rIP.txt")
# oneResult50("n50BW50\\rDSIS_MII.txt")
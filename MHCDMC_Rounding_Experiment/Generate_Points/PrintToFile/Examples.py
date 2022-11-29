import os.path

from generate_points.GenerateLocation import generateUserLocationRandomly
from generate_points.GenerateLocation import generateServerLocationByNeibor
from generate_points.GenerateRequestAndCapacity import generateUserRequest
from generate_points.GenerateRequestAndCapacity import generateServerCapacity
from draw_points.DrawPoints import draw_all_points
from PrintToFile.PrintDataToFile import print_data_to_file

def one_instance():
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




one_instance()


